"""
环节1：分子标准化流程

严格按设计文档"1.2 标准化流程"的顺序执行，顺序不能乱：
  1. RDKit MolFromSmiles(sanitize=True)，失败的进隔离表，不静默丢弃
  2. 去盐 + 取最大有机片段
  3. 电荷中和（不在这一步固定质子化态，只中和形式电荷）
  4. 互变异构体规范化
  5. 保留所有立体中心
  6. 生成 InChIKey14（去立体，骨架聚合用）和完整 InChIKey（精确去重用）两套主键

用法：
    from core.standardize import standardize_batch
    result = standardize_batch(["CCO", "not-a-smiles", "c1ccccc1.[Na+]"])
    result.curated     # -> DataFrame，标准化成功的行
    result.quarantine  # -> DataFrame，失败的行 + 失败原因，供人工复核
"""
from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import inchi
from rdkit.Chem.MolStandardize import rdMolStandardize

# RDKit 默认会把每个 sanitize 警告都打到 stderr，批量跑几千个分子时输出会被刷屏。
# 失败/警告信息我们自己在 StandardizeResult.quarantine 里结构化记录，不需要 RDKit 自己再打一遍。
RDLogger.DisableLog("rdApp.*")

_FRAGMENT_CHOOSER = rdMolStandardize.LargestFragmentChooser()
_UNCHARGER = rdMolStandardize.Uncharger()
_TAUTOMER_ENUMERATOR = rdMolStandardize.TautomerEnumerator()


@dataclass
class StandardizeResult:
    curated: pd.DataFrame = field(default_factory=pd.DataFrame)
    quarantine: pd.DataFrame = field(default_factory=pd.DataFrame)

    @property
    def success_rate(self) -> float:
        total = len(self.curated) + len(self.quarantine)
        return len(self.curated) / total if total else 0.0


def standardize_one(raw_smiles: str) -> dict:
    """
    对单个 SMILES 跑完整标准化流程。
    返回一个 dict：成功时包含 standardized_smiles/inchikey/inchikey14 等字段且 ok=True；
    失败时 ok=False 并带 fail_stage/fail_reason，调用方据此决定是否进隔离表。
    """
    record: dict = {"raw_smiles": raw_smiles, "ok": False}

    # 1. 解析 + sanitize
    try:
        mol = Chem.MolFromSmiles(raw_smiles, sanitize=True)
    except Exception as exc:  # RDKit 某些异常是 C++ 抛出的非标准异常，用宽 except 兜底
        mol = None
        record["fail_stage"] = "parse"
        record["fail_reason"] = str(exc)
    if mol is None:
        record.setdefault("fail_stage", "parse")
        record.setdefault("fail_reason", "MolFromSmiles returned None (invalid SMILES or sanitize failed)")
        return record

    # 2. 去盐 + 取最大有机片段
    try:
        mol = _FRAGMENT_CHOOSER.choose(mol)
    except Exception as exc:
        record["fail_stage"] = "largest_fragment"
        record["fail_reason"] = str(exc)
        return record

    # 3. 电荷中和（不固定质子化态，只是把能中和的形式电荷中和掉）
    try:
        mol = _UNCHARGER.uncharge(mol)
    except Exception as exc:
        record["fail_stage"] = "uncharge"
        record["fail_reason"] = str(exc)
        return record

    # 4. 互变异构体规范化——这步比去盐重要，吡唑/酰胺互变会让同一化合物出现两个不同InChIKey
    try:
        mol = _TAUTOMER_ENUMERATOR.Canonicalize(mol)
    except Exception as exc:
        record["fail_stage"] = "tautomer_canonicalize"
        record["fail_reason"] = str(exc)
        return record

    # 5. 立体化学：这里全程没有调用任何会清除立体标记的函数(如 RemoveStereochemistry)，
    #    上面三步(LargestFragmentChooser/Uncharger/TautomerEnumerator)都是保留立体中心的操作。
    #    显式做一次合法性检查，确认手性标记在标准化后仍然有效。
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

    # 6. 双主键：InChIKey(全串，精确去重) + InChIKey14(去立体，骨架聚合)
    try:
        full_inchi = inchi.MolToInchi(mol)
        full_key = inchi.InchiToInchiKey(full_inchi)
        key14 = full_key.split("-")[0] if full_key else None
    except Exception as exc:
        record["fail_stage"] = "inchi"
        record["fail_reason"] = str(exc)
        return record

    canonical_smiles = Chem.MolToSmiles(mol)
    record.update(
        {
            "ok": True,
            "standardized_smiles": canonical_smiles,
            "inchikey": full_key,
            "inchikey14": key14,
            "num_stereo_centers": len(Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)),
            "formal_charge": Chem.GetFormalCharge(mol),
        }
    )
    return record


def standardize_batch(raw_smiles_list: list[str]) -> StandardizeResult:
    """批量标准化，自动把失败的分子分流进 quarantine，不静默丢弃（环节1.2 第1步的要求）。"""
    rows = [standardize_one(smi) for smi in raw_smiles_list]
    df = pd.DataFrame(rows)
    if df.empty:
        return StandardizeResult()
    # 全部成功（或全部失败）时，df 里可能根本不存在 fail_stage/fail_reason 列——
    # list[dict]里没有任何一行带这个key，pandas就不会生成该列，下面按列名取子集会抛KeyError。
    # 这个bug是跑workflows/run_discovery_round.py端到端集成时才暴露出来的，单测没覆盖到。
    for col in ("fail_stage", "fail_reason"):
        if col not in df.columns:
            df[col] = None
    curated = df[df["ok"]].drop(columns=["ok", "fail_stage", "fail_reason"], errors="ignore").reset_index(drop=True)
    quarantine = df[~df["ok"]][["raw_smiles", "fail_stage", "fail_reason"]].reset_index(drop=True)
    return StandardizeResult(curated=curated, quarantine=quarantine)


if __name__ == "__main__":
    # 冒烟测试：包含一个正常分子、一个带盐的分子、一个非法SMILES
    demo = [
        "COc1cc(N(C)CCN(C)C)c(NC(=O)C=C)cc1Nc1nccc(-c2cn(C)c3ccccc23)n1",  # 奥希替尼(已验证分子式正确)
        "c1ccccc1.[Na+].[Cl-]",  # 苯 + 氯化钠盐，应该只剩苯
        "this-is-not-a-smiles(((",
    ]
    result = standardize_batch(demo)
    print("=== curated ===")
    print(result.curated.to_string())
    print("\n=== quarantine ===")
    print(result.quarantine.to_string())
    print(f"\nsuccess_rate = {result.success_rate:.2%}")
