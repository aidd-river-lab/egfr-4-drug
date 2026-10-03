"""
core/ternary_complex.py 的回归测试，重点是 evaluate_hook_effect_risk() ——
双功能降解剂公认的真实药理学现象，不需要三元复合物3D预测就能验证这条纯数学判据。
"""
from core.ternary_complex import TernaryComplexResult, evaluate_hook_effect_risk, run_ternary_complex_stub


def test_hook_effect_not_triggered_within_safe_dose_range():
    result = evaluate_hook_effect_risk(warhead_kd_nm=50, e3_ligand_kd_nm=20, dose_range_nm=(10, 100))
    assert result["hook_effect_risk"] is False


def test_hook_effect_triggered_when_dose_exceeds_threshold():
    # 阈值 = 10 * max(50, 20) = 500；剂量上限1000超过阈值
    result = evaluate_hook_effect_risk(warhead_kd_nm=50, e3_ligand_kd_nm=20, dose_range_nm=(100, 1000))
    assert result["hook_effect_risk"] is True


def test_hook_effect_threshold_uses_the_larger_kd():
    """阈值取两个Kd里更大的那个(limiting_kd)，不是取平均或更小值——这是药理学上的保守做法。"""
    result = evaluate_hook_effect_risk(warhead_kd_nm=200, e3_ligand_kd_nm=20, dose_range_nm=(10, 100))
    assert result["limiting_kd_nm"] == 200
    assert result["hook_risk_threshold_nm"] == 2000


def test_hook_effect_boundary_is_inclusive():
    # 剂量上限恰好等于阈值，算作进入风险区间(>=，不是>)
    result = evaluate_hook_effect_risk(warhead_kd_nm=50, e3_ligand_kd_nm=20, dose_range_nm=(10, 500))
    assert result["hook_effect_risk"] is True


def test_safe_margin_fold_is_threshold_over_dose_low():
    result = evaluate_hook_effect_risk(warhead_kd_nm=50, e3_ligand_kd_nm=20, dose_range_nm=(10, 100))
    assert result["safe_margin_fold"] == 50.0  # 500 / 10


def test_run_ternary_complex_stub_is_honest_about_not_running():
    result = run_ternary_complex_stub("demo_warhead", "lenalidomide_n_linked", "peg2")
    assert result["ok"] is False
    assert "PRosettaC" in result["reason"] or "Rosetta" in result["reason"]


def test_ternary_complex_result_leaves_unmeasured_fields_none():
    """真实三元复合物数据(cooperativity_alpha/DC50/Dmax)在没有3D预测工具时必须是None，不能有默认值填充。"""
    result = TernaryComplexResult(
        warhead_id="w", e3_ligand_id="e", linker_id="l", warhead_kd_nm=50.0, e3_ligand_kd_nm=20.0
    )
    assert result.cooperativity_alpha is None
    assert result.predicted_dc50_nm is None
    assert result.predicted_dmax_pct is None
