"""
core/selectivity.py 的符号约定回归测试。

背景：第一版docstring写反过一次方向(误写成"ΔΔΔG>0表示对突变体更有利")，
是在验证demo输出之后才发现并改正的。这里用无歧义的合成场景把方向钉死，
防止以后谁改代码时又把符号绕反。
"""
import math

from core.selectivity import (
    compute_dddg_from_fep_maps,
    dddg_to_fold_selectivity,
    propagate_subtraction_uncertainty,
)
from core.fep import FEPMapResult, FEPPerturbation


def test_negative_dddg_means_mutant_favored():
    """ΔΔΔG < 0 必须对应 fold_selectivity > 1(更偏向突变体)。"""
    fold = dddg_to_fold_selectivity(-3.0)
    assert fold > 1.0


def test_positive_dddg_means_wt_favored():
    """ΔΔΔG > 0 必须对应 fold_selectivity < 1(更偏向WT)。"""
    fold = dddg_to_fold_selectivity(3.0)
    assert fold < 1.0


def test_zero_dddg_means_no_selectivity():
    fold = dddg_to_fold_selectivity(0.0)
    assert math.isclose(fold, 1.0, rel_tol=1e-9)


def test_unambiguous_synthetic_scenario_mutant_strongly_favored():
    """
    构造一个方向无歧义的场景：同一个改造(L1->L2)在突变体上结合大幅变好(RBFE=-3.0)，
    在WT上完全不变(RBFE=0.0)。这种改造显然应该让L2对突变体更有选择性，fold必须>>1。
    """
    mutant_map = FEPMapResult(
        genotype="del19_C797S",
        perturbations=[FEPPerturbation("L1", "L2", predicted_ddg_kcal_mol=-3.0, uncertainty_kcal_mol=0.3)],
    )
    wt_map = FEPMapResult(
        genotype="EGFR_WT",
        perturbations=[FEPPerturbation("L1", "L2", predicted_ddg_kcal_mol=0.0, uncertainty_kcal_mol=0.3)],
    )
    result = compute_dddg_from_fep_maps("L1", "L2", "del19_C797S", mutant_map, wt_map)
    assert result["ok"] is True
    assert result["dddg_kcal_mol"] < 0
    assert result["fold_selectivity_pred"] > 100  # exp(3.0 / (R*T)) ~= 158


def test_missing_perturbation_pair_reported_honestly():
    mutant_map = FEPMapResult(genotype="del19_C797S", perturbations=[])
    wt_map = FEPMapResult(genotype="EGFR_WT", perturbations=[])
    result = compute_dddg_from_fep_maps("L1", "L2", "del19_C797S", mutant_map, wt_map)
    assert result["ok"] is False


def test_uncertainty_propagation_is_quadrature_sum():
    u = propagate_subtraction_uncertainty(0.3, 0.4)
    assert math.isclose(u, 0.5, rel_tol=1e-9)  # 3-4-5 直角三角形，凑整好验证
