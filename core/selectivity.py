"""
环节5：选择性引擎——设计文档原话"整条管线最该投入算力的地方"

SelectivityResult 的字段完全照抄设计文档"5.4 代码落地"给出的接口形状，不要改字段名，
下游(core/decide.py 的多目标打分)按这个结构读取。

本模块真实、可测试的部分是热力学转换和误差传播：
  - ΔΔΔG(kcal/mol) -> 倍数选择性的转换(玻尔兹曼关系，温度298K)
  - 两次RBFE相减时的不确定度传播(标准误差传播公式)
这两步都是纯数学，不需要真的跑过FEP就能验证对不对。
需要真实FEP/MD计算的部分委托给 core/fep.py（已经在那边诚实标注为未执行）。
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Literal

from core.fep import FEPMapResult, FEPPerturbation

R_KCAL_PER_MOL_K = 1.987204e-3  # 气体常数，kcal/(mol·K)
TEMPERATURE_K = 298.15


@dataclass
class SelectivityResult:
    """字段形状完全照抄设计文档环节5.4，不要改名字。"""

    ligand_id: str
    ddg_target: dict[str, float] = field(default_factory=dict)  # 每个目标基因型的 ΔG
    ddg_wt: float | None = None
    dddg: dict[str, float] = field(default_factory=dict)  # 选择性差值，排序用主指标
    fold_selectivity_pred: dict[str, float] = field(default_factory=dict)
    uncertainty: dict[str, float] = field(default_factory=dict)
    method: Literal["fep", "mmgbsa", "ml", "docking"] = "fep"
    structures_used: list[str] = field(default_factory=list)


def dddg_to_fold_selectivity(dddg_kcal_mol: float, temperature_k: float = TEMPERATURE_K) -> float:
    """
    ΔΔΔG(kcal/mol) 转换成倍数选择性，基于 ΔG = -RT ln(K) 的玻尔兹曼关系。

    符号约定(容易搞反，务必对照这段再用)：RBFE(L1->L2) = ΔG_L2 - ΔG_L1，按FEP惯例
    "更负 = 结合更有利"。ΔΔΔG = RBFE_突变体 - RBFE_WT。
    **ΔΔΔG < 0** 表示"L1->L2这个改造带来的结合增益，在突变体里比在WT里更大"，
    即 L2 相对 L1 对突变体更有选择性偏好，fold > 1。
    ΔΔΔG > 0 则相反：这个改造让分子更偏向WT，选择性变差，fold < 1。
    已经用一个无歧义场景验证过这个方向是对的(见 tests/test_selectivity.py)。
    """
    return math.exp(-dddg_kcal_mol / (R_KCAL_PER_MOL_K * temperature_k))


def propagate_subtraction_uncertainty(u1: float, u2: float) -> float:
    """两个独立估计值相减时的标准误差传播：sigma = sqrt(sigma1^2 + sigma2^2)。"""
    return math.sqrt(u1**2 + u2**2)


def compute_dddg_from_fep_maps(
    ligand_a: str,
    ligand_b: str,
    genotype: str,
    mutant_map: FEPMapResult,
    wt_map: FEPMapResult,
) -> dict:
    """
    环节5.1的炼金术双炼金术：在突变体和WT上跑同一组配体微扰(ligand_a -> ligand_b)，相减。
    要求两张FEP map里都能找到完全一样的(ligand_a, ligand_b)微扰对，否则无法相减——
    这是方法学上的硬要求，不是本实现偷懒。
    """

    def _find(fep_map: FEPMapResult):
        for p in fep_map.perturbations:
            if p.ligand_a == ligand_a and p.ligand_b == ligand_b:
                return p
        return None

    mut_pert = _find(mutant_map)
    wt_pert = _find(wt_map)
    if mut_pert is None or wt_pert is None:
        return {
            "ok": False,
            "reason": f"在突变体({'找到' if mut_pert else '没找到'})和WT({'找到' if wt_pert else '没找到'})"
            f"map里必须都有完全一样的微扰对({ligand_a}->{ligand_b})才能算ΔΔΔG",
        }

    dddg = mut_pert.predicted_ddg_kcal_mol - wt_pert.predicted_ddg_kcal_mol
    uncertainty = propagate_subtraction_uncertainty(mut_pert.uncertainty_kcal_mol, wt_pert.uncertainty_kcal_mol)
    fold = dddg_to_fold_selectivity(dddg)

    return {
        "ok": True,
        "genotype": genotype,
        "ligand_pair": (ligand_a, ligand_b),
        "dddg_kcal_mol": round(dddg, 3),
        "uncertainty_kcal_mol": round(uncertainty, 3),
        "fold_selectivity_pred": round(fold, 1),
        "interpretation": f"{ligand_b} 相对 {ligand_a} 在 {genotype} 上比在WT上多出约 {fold:.1f} 倍的相对亲和力优势",
    }


if __name__ == "__main__":
    # 合成FEP数据演示(不是真实计算结果，用来验证ΔΔΔG计算+热力学转换+误差传播逻辑)
    mutant_map = FEPMapResult(
        genotype="del19_C797S",
        perturbations=[FEPPerturbation("compound_A", "compound_B", predicted_ddg_kcal_mol=-1.5, uncertainty_kcal_mol=0.3)],
    )
    wt_map = FEPMapResult(
        genotype="EGFR_WT",
        perturbations=[FEPPerturbation("compound_A", "compound_B", predicted_ddg_kcal_mol=0.2, uncertainty_kcal_mol=0.25)],
    )

    result = compute_dddg_from_fep_maps("compound_A", "compound_B", "del19_C797S", mutant_map, wt_map)
    print("=== ΔΔΔG 计算结果 ===")
    print(result)

    print("\n=== 构造一个完整的 SelectivityResult ===")
    sel = SelectivityResult(
        ligand_id="compound_B",
        ddg_target={"del19_C797S": -1.5},
        ddg_wt=0.2,
        dddg={"del19_C797S": result["dddg_kcal_mol"]},
        fold_selectivity_pred={"del19_C797S": result["fold_selectivity_pred"]},
        uncertainty={"del19_C797S": result["uncertainty_kcal_mol"]},
        method="fep",
        structures_used=["del19_C797S_cluster1", "EGFR_WT_cluster1"],
    )
    print(sel)
