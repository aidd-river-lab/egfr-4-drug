"""
环节3：类似物枚举（按口袋区域定向的 R 基团组合，不是随机加官能团）

流程（环节3.2）：
  1. 从 config/rgroup_libraries/scaffolds.yaml 读取骨架模板和可变位点
  2. 每个位点按 region 去 config/rgroup_libraries/pocket_regions.yaml 找对应片段库
  3. 组合枚举（用 RDKit molzip 按同位素标记的哑原子做片段连接，不是字符串拼接）
  4. 立即跑合成可及性粗过滤（本模块只接口对接 core/synthesis.py，不做过滤本身）

用法：
    from core.enumerate import enumerate_from_scaffold
    df = enumerate_from_scaffold("demo_aminopyrimidine_biphenyl")
    df.columns -> ['smiles', 'scaffold_id', 'R1_name', 'R2_name', 'R3_name']
"""
from __future__ import annotations

import itertools
from pathlib import Path

import pandas as pd
import yaml
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

CONFIG_DIR = Path(__file__).resolve().parent.parent / "config" / "rgroup_libraries"


def _load_yaml(name: str) -> dict:
    with open(CONFIG_DIR / name, encoding="utf-8") as f:
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


def load_scaffold(scaffold_id: str) -> dict:
    """从 scaffolds.yaml 里取出指定 id 的骨架定义。"""
    data = _load_yaml("scaffolds.yaml")
    for tpl in data["templates"]:
        if tpl["id"] == scaffold_id:
            return tpl
    available = [t["id"] for t in data["templates"]]
    raise KeyError(f"未找到骨架 {scaffold_id!r}，可用: {available}")


def load_region_fragments(region: str) -> list[dict]:
    """从 pocket_regions.yaml 里取出指定口袋区域的全部片段。"""
    data = _load_yaml("pocket_regions.yaml")
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


def enumerate_from_scaffold(scaffold_id: str, max_combinations: int | None = None) -> pd.DataFrame:
    """
    对指定骨架做全组合枚举（笛卡尔积）。
    小规模的 pocket_regions 库（几到十几个/区域）做全枚举是可行的；
    如果以后接入商业R基团库（成百上千/区域），要用 max_combinations 做随机抽样封顶，
    或者换成环节3.3 的"受约束生成"而不是穷举。
    """
    tpl = load_scaffold(scaffold_id)
    attachment_points = tpl["attachment_points"]

    # 每个位点的候选片段列表：[(isotope, fragment_dict), ...] 的笛卡尔积
    per_point_fragments = []
    for pt in attachment_points:
        frags = load_region_fragments(pt["region"])
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
