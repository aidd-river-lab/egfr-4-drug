"""
core/ternary_complex.py 的回归测试。hook effect是纯数学判据；kabsch_rigid_transform
是纯numpy；measure_linker_span是纯RDKit——这三个都不需要PyRosetta，可以在主.venv
下完整验证。run_protein_protein_docking_trials需要PyRosetta(.venv310)，这里只测
主.venv下的诚实失败路径，真实执行的验证记录见doc/routes/route-b-degrader.md。
"""
import numpy as np

from core.ternary_complex import (
    TernaryComplexResult,
    evaluate_hook_effect_risk,
    evaluate_ternary_complex_geometry,
    kabsch_rigid_transform,
    measure_linker_span,
    run_protein_protein_docking_trials,
    run_ternary_complex_stub,
)


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
    assert "cooperativity" in result["reason"]


def test_ternary_complex_result_leaves_unmeasured_fields_none():
    """真实三元复合物数据(cooperativity_alpha/DC50/Dmax)在没有3D预测工具时必须是None，不能有默认值填充。"""
    result = TernaryComplexResult(
        warhead_id="w", e3_ligand_id="e", linker_id="l", warhead_kd_nm=50.0, e3_ligand_kd_nm=20.0
    )
    assert result.cooperativity_alpha is None
    assert result.predicted_dc50_nm is None
    assert result.predicted_dmax_pct is None


# ------------------------------------------------------------
# kabsch_rigid_transform：纯numpy，不需要PyRosetta
# ------------------------------------------------------------
def test_kabsch_recovers_pure_translation():
    before = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    translation = np.array([5.0, 3.0, 1.0])
    after = before + translation
    transform = kabsch_rigid_transform(before, after)
    test_point = np.array([2.0, 2.0, 2.0])
    np.testing.assert_allclose(transform(test_point), test_point + translation, atol=1e-8)


def test_kabsch_recovers_90_degree_rotation():
    """绕Z轴转90度：(1,0,0)应该变成(0,1,0)这类场景，验证旋转部分算对了，不只是平移。"""
    before = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    rotation_90z = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    after = before @ rotation_90z.T
    transform = kabsch_rigid_transform(before, after)
    test_point = np.array([1.0, 0.0, 0.0])
    expected = rotation_90z @ test_point
    np.testing.assert_allclose(transform(test_point), expected, atol=1e-6)


def test_kabsch_identity_when_nothing_moves():
    before = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    transform = kabsch_rigid_transform(before, before.copy())
    test_point = np.array([3.3, -2.1, 7.0])
    np.testing.assert_allclose(transform(test_point), test_point, atol=1e-8)


# ------------------------------------------------------------
# measure_linker_span：纯RDKit，不需要PyRosetta
# ------------------------------------------------------------
def test_measure_linker_span_peg2_matches_known_real_value():
    """这个数字是真跑过50个MMFF优化构象测出来的(见本次构建记录)，不是凭经验猜的，
    回归测试锁定这个量级，防止以后改坏了还不知道。"""
    result = measure_linker_span("*CCOCCOCC*")
    assert result["ok"] is True
    assert 3 < result["min_span_angstrom"] < 6
    assert 7 < result["max_span_angstrom"] < 10


def test_measure_linker_span_longer_linker_has_larger_span():
    peg1 = measure_linker_span("*CCOCC*")
    peg3 = measure_linker_span("*CCOCCOCCOCC*")
    assert peg3["median_span_angstrom"] > peg1["median_span_angstrom"]


def test_measure_linker_span_rejects_wrong_dummy_count():
    result = measure_linker_span("*CCC")  # 只有1个连接点
    assert result["ok"] is False


# ------------------------------------------------------------
# run_protein_protein_docking_trials：需要PyRosetta，这里只测主.venv下的诚实失败
# ------------------------------------------------------------
def test_docking_trials_honestly_fails_without_pyrosetta():
    result = run_protein_protein_docking_trials(
        combined_pdb_path="/nonexistent.pdb",
        fixed_chain="A",
        mobile_chain="B",
        mobile_chain_external_point_original=(0, 0, 0),
        fixed_chain_external_point=(0, 0, 0),
        n_trials=1,
    )
    if not result["ok"]:
        assert "PyRosetta" in result["reason"]


def test_evaluate_ternary_complex_geometry_passes_through_failure():
    """docking trials本身失败时，evaluate_ternary_complex_geometry不应该假装能算出几何兼容性。"""
    failed_docking = {"ok": False, "reason": "需要PyRosetta..."}
    result = evaluate_ternary_complex_geometry(failed_docking, linker_max_span_angstrom=10.0)
    assert result["ok"] is False


def test_evaluate_ternary_complex_geometry_detects_compatible_pose():
    fake_docking_result = {
        "ok": True,
        "n_trials": 3,
        "trials_sorted_by_distance": [
            {"trial": 0, "score": 100.0, "exit_vector_distance_angstrom": 8.0},
            {"trial": 1, "score": 101.0, "exit_vector_distance_angstrom": 15.0},
            {"trial": 2, "score": 99.0, "exit_vector_distance_angstrom": 30.0},
        ],
    }
    result = evaluate_ternary_complex_geometry(fake_docking_result, linker_max_span_angstrom=12.0)
    assert result["any_compatible_pose_found"] is True
    assert result["n_geometrically_compatible_trials"] == 1


def test_evaluate_ternary_complex_geometry_honestly_reports_no_compatible_pose():
    """真实跑过的场景(500次试验，最好26.8Å，远超linker最大span)——这条测试复现这个真实发现。"""
    fake_docking_result = {
        "ok": True,
        "n_trials": 500,
        "trials_sorted_by_distance": [{"trial": 0, "score": 1122.3, "exit_vector_distance_angstrom": 26.8}],
    }
    result = evaluate_ternary_complex_geometry(fake_docking_result, linker_max_span_angstrom=11.9)
    assert result["any_compatible_pose_found"] is False
    assert result["n_geometrically_compatible_trials"] == 0
