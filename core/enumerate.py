"""
环节3：类似物枚举

两种拓扑：
  1. 固定骨架 + R基团(路线A/C用)：从 scaffolds.yaml 读取骨架模板和可变位点，
     每个位点按 region 去 pocket_regions.yaml 找对应片段库，组合枚举。
     用法: enumerate_from_scaffold("demo_aminopyrimidine_biphenyl")
     df.columns -> ['smiles', 'scaffold_id', 'R1_name', 'R2_name', 'R3_name']

  2. warhead-linker-E3三组分(路线B双功能降解剂专用)：没有"谁是骨架"的区分，
     三个独立片段依次相连。用法: enumerate_bifunctional_library(config_dir=...)
     df.columns -> ['smiles', 'warhead_name', 'linker_name', 'e3_ligand_name']

两种拓扑都用 RDKit molzip 按同位素标记的哑原子做片段连接(不是字符串拼接)，
都做了同一类连通性校验(见 assemble()/enumerate_bifunctional() 内部注释)。
"""
from __future__ import annotations

import itertools
from pathlib import Path

import pandas as pd
import yaml
from rdkit import Chem, RDLogger

from core.route_config import route_config_dir

RDLogger.DisableLog("rdApp.*")


def _load_yaml(name: str, config_dir: Path | None = None) -> dict:
    config_dir = config_dir or route_config_dir()
    with open(config_dir / "rgroup_libraries" / name, encoding="utf-8") as f:
        return yaml.safe_load(f)


def _mol_with_isotope_dummy(smiles: str, force_isotope: int | None = None) -> Chem.Mol:
    """
    解析片段 SMILES。pocket_regions.yaml 里的片段统一用裸 "*" 标记连接点（同一份库要
    被复用到不同的骨架位点上，不能在配置里写死某个isotope），这里按调用方传入的
    force_isotope 把哑原子的同位素标记设成骨架要求的编号，molzip 才能正确配对——
    如果不做这一步，molzip 会因为找不到匹配的同位素标记而把片段原样留成游离的"."分子，
    SanitizeMol 不会报错，只会生成一个看起来正常实际上根本没拼起来的多组分SMILES。
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"无法解析片段 SMILES: {smiles!r}")
    if force_isotope is not None:
        dummy_atoms = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
        if len(dummy_atoms) != 1:
            raise ValueError(f"片段 {smiles!r} 应该恰好有1个连接点(*)，实际有{len(dummy_atoms)}个")
        dummy_atoms[0].SetIsotope(force_isotope)
    return mol


def _mol_with_isotope_dummies(smiles: str, isotopes: list[int]) -> Chem.Mol:
    """
    _mol_with_isotope_dummy() 的多连接点版本，linker片段(两端各一个连接点)专用。
    不复用同一个函数名/放宽原来"恰好1个"的校验，是为了不削弱骨架+R基团路径上
    那条已经抓到过真实bug的连通性校验——这里是另一类拓扑(linker两端分别接不同的
    邻居)，该严格的地方不能因为复用代码而放松。
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"无法解析片段 SMILES: {smiles!r}")
    dummy_atoms = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
    if len(dummy_atoms) != len(isotopes):
        raise ValueError(f"片段 {smiles!r} 应该有{len(isotopes)}个连接点(*)，实际有{len(dummy_atoms)}个")
    for atom, isotope in zip(dummy_atoms, isotopes):
        atom.SetIsotope(isotope)
    return mol


def enumerate_bifunctional(warhead_smiles: str, linker_smiles: str, e3_ligand_smiles: str) -> str | None:
    """
    路线B专用：warhead(1个连接点) - linker(2个连接点) - e3_ligand(1个连接点)三组分拼接。
    和 assemble() 的固定骨架+R基团拓扑不同，这里三个片段地位对等，没有谁是"骨架"；
    linker的两端分别对接warhead和e3_ligand。复用同一套molzip+连通性校验逻辑——
    "哑原子同位素不匹配导致molzip静默留下断裂分子"这个风险在三组分拓扑下同样存在，
    甚至更容易出错(两个独立的同位素配对，而不是骨架那种"学生"式的单侧配对)。
    返回标准化SMILES；拼接/sanitize失败返回None（不抛异常中断整批枚举）。
    """
    try:
        w = _mol_with_isotope_dummy(warhead_smiles, force_isotope=1)
        link = _mol_with_isotope_dummies(linker_smiles, isotopes=[1, 2])
        e3 = _mol_with_isotope_dummy(e3_ligand_smiles, force_isotope=2)
        combined = Chem.CombineMols(Chem.CombineMols(w, link), e3)
        params = Chem.MolzipParams()
        params.label = Chem.MolzipLabel.Isotope
        result = Chem.molzip(combined, params)
        Chem.SanitizeMol(result)
        if len(Chem.GetMolFrags(result)) != 1:
            return None
        if any(a.GetAtomicNum() == 0 for a in result.GetAtoms()):
            return None
        return Chem.MolToSmiles(result)
    except Exception:
        return None


def load_bifunctional_fragments(config_dir: Path | None = None) -> dict:
    """从 warhead_library.yaml / linkers.yaml / e3_ligands.yaml 读取路线B的三类片段库。"""
    warheads = _load_yaml("warhead_library.yaml", config_dir)["warheads"]["fragments"]
    linkers = _load_yaml("linkers.yaml", config_dir)["linkers"]["fragments"]
    e3_ligands = _load_yaml("e3_ligands.yaml", config_dir)["e3_ligands"]["fragments"]
    return {"warheads": warheads, "linkers": linkers, "e3_ligands": e3_ligands}


def enumerate_bifunctional_library(config_dir: Path | None = None, max_combinations: int | None = None) -> pd.DataFrame:
    """对 warhead x linker x e3_ligand 做全组合枚举，返回形状类似 enumerate_from_scaffold() 的DataFrame。"""
    config_dir = config_dir or route_config_dir()
    libs = load_bifunctional_fragments(config_dir)

    rows = []
    combo_iter = itertools.product(libs["warheads"], libs["linkers"], libs["e3_ligands"])
    for i, (w, link, e3) in enumerate(combo_iter):
        if max_combinations is not None and i >= max_combinations:
            break
        smiles = enumerate_bifunctional(w["smiles"], link["smiles"], e3["smiles"])
        if smiles is None:
            continue
        rows.append(
            {"smiles": smiles, "warhead_name": w["name"], "linker_name": link["name"], "e3_ligand_name": e3["name"]}
        )

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.drop_duplicates(subset="smiles").reset_index(drop=True)
    return df


def load_scaffold(scaffold_id: str, config_dir: Path | None = None) -> dict:
    """从 scaffolds.yaml 里取出指定 id 的骨架定义。"""
    data = _load_yaml("scaffolds.yaml", config_dir)
    for tpl in data["templates"]:
        if tpl["id"] == scaffold_id:
            return tpl
    available = [t["id"] for t in data["templates"]]
    raise KeyError(f"未找到骨架 {scaffold_id!r}，可用: {available}")


def load_region_fragments(region: str, config_dir: Path | None = None) -> list[dict]:
    """从 pocket_regions.yaml 里取出指定口袋区域的全部片段。"""
    data = _load_yaml("pocket_regions.yaml", config_dir)
    if region not in data:
        raise KeyError(f"pocket_regions.yaml 里没有区域 {region!r}，可用: {list(data.keys())}")
    return data[region]["fragments"]


def assemble(core_smiles: str, fragment_smiles_by_isotope: dict[int, str]) -> str | None:
    """
    用 RDKit molzip 把骨架和一组按同位素标记匹配的片段拼成最终分子。
    fragment_smiles_by_isotope: {1: "[1*]S(C)(=O)=O", 2: "[2*]C(F)(F)F", ...}
    返回标准化后的 SMILES；拼接/sanitize 失败返回 None（不抛异常中断整批枚举）。
    """
    try:
        combined = _mol_with_isotope_dummy(core_smiles)  # 骨架本身的哑原子已经带正确编号，不需要force
        for isotope, smi in fragment_smiles_by_isotope.items():
            frag = _mol_with_isotope_dummy(smi, force_isotope=isotope)
            combined = Chem.CombineMols(combined, frag)
        params = Chem.MolzipParams()
        params.label = Chem.MolzipLabel.Isotope
        result = Chem.molzip(combined, params)
        Chem.SanitizeMol(result)
        if len(Chem.GetMolFrags(result)) != 1:
            # molzip 在找不到匹配同位素标记时不会报错，只会把哑原子悄悄留在原地，
            # 产出一个看似合法、实际上骨架和片段根本没连起来的多组分分子（本模块开发时
            # 真实踩过这个坑）。这里做硬校验，宁可整条丢弃也不能让断裂分子混进枚举结果。
            return None
        if any(a.GetAtomicNum() == 0 for a in result.GetAtoms()):
            return None  # 还有残留哑原子，说明某个连接点没有被正确消耗
        return Chem.MolToSmiles(result)
    except Exception:
        return None


def enumerate_from_scaffold(
    scaffold_id: str, max_combinations: int | None = None, config_dir: Path | None = None
) -> pd.DataFrame:
    """
    对指定骨架做全组合枚举（笛卡尔积）。
    小规模的 pocket_regions 库（几到十几个/区域）做全枚举是可行的；
    如果以后接入商业R基团库（成百上千/区域），要用 max_combinations 做随机抽样封顶，
    或者换成环节3.3 的"受约束生成"而不是穷举。

    config_dir: 不传则默认路线C的配置；路线A/B传各自的 route_config_dir(...)。
    """
    tpl = load_scaffold(scaffold_id, config_dir)
    attachment_points = tpl["attachment_points"]

    # 每个位点的候选片段列表：[(isotope, fragment_dict), ...] 的笛卡尔积
    per_point_fragments = []
    for pt in attachment_points:
        frags = load_region_fragments(pt["region"], config_dir)
        per_point_fragments.append([(pt["isotope"], pt["label"], f) for f in frags])

    rows = []
    combo_iter = itertools.product(*per_point_fragments)
    for i, combo in enumerate(combo_iter):
        if max_combinations is not None and i >= max_combinations:
            break
        frag_by_isotope = {isotope: f["smiles"] for isotope, _label, f in combo}
        smiles = assemble(tpl["core_smiles"], frag_by_isotope)
        if smiles is None:
            continue
        row = {"smiles": smiles, "scaffold_id": scaffold_id}
        for isotope, label, f in combo:
            row[f"{label}_name"] = f["name"]
        rows.append(row)

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.drop_duplicates(subset="smiles").reset_index(drop=True)
    return df


if __name__ == "__main__":
    df = enumerate_from_scaffold("demo_aminopyrimidine_biphenyl")
    print(f"枚举出 {len(df)} 个去重后的候选分子（骨架: demo_aminopyrimidine_biphenyl）")
    print(df.head(10).to_string())

    print("\n=== 路线B：warhead-linker-E3三组分拼接演示 ===")
    bifunc_df = enumerate_bifunctional_library(config_dir=route_config_dir("route_b_degrader"))
    print(f"枚举出 {len(bifunc_df)} 个双功能降解剂候选")
    print(bifunc_df.head(5).to_string())
