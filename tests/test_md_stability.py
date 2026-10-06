"""
core/md_stability.py 的回归测试。evaluate_pass_criteria()之前只在__main__里用
合成数据演示过，没有真正进pytest——这里把那个演示补成真正的回归测试。
run_protein_ligand_complex_md()需要openmm+openff-toolkit+openmmforcefields+
真实AmberTools(装在一个独立的conda环境，不是.venv/.venv310，见该函数docstring)，
这里只测主.venv下的诚实失败路径；真实执行结果记录见doc/04-funnel-docking.md。
"""
from core.md_stability import MDStabilityResult, evaluate_pass_criteria, run_md_stability_stub, run_protein_ligand_complex_md


def test_run_md_stability_stub_is_honest_about_not_running():
    result = run_md_stability_stub("DEMO-001", "receptor_clean.pdb")
    assert result["ok"] is False
    assert "GPU" in result["reason"]


def test_evaluate_pass_criteria_good_result_passes():
    good = MDStabilityResult(
        "DEMO-GOOD", ligand_rmsd_last_15ns_mean_angstrom=1.2, hinge_hbond_occupancy_pct=85,
        target_anchor_occupancy_pct=55, md_length_ns=20, n_replicas=3,
    )
    result = evaluate_pass_criteria(good)
    assert result["rmsd_ok"] is True
    assert result["hinge_hbond_ok"] is True
    assert result["target_anchor_ok"] is True
    assert result["overall_pass"] is True


def test_evaluate_pass_criteria_bad_result_fails_all_three():
    bad = MDStabilityResult(
        "DEMO-BAD", ligand_rmsd_last_15ns_mean_angstrom=4.1, hinge_hbond_occupancy_pct=30,
        target_anchor_occupancy_pct=10, md_length_ns=20, n_replicas=3,
    )
    result = evaluate_pass_criteria(bad)
    assert result["rmsd_ok"] is False
    assert result["hinge_hbond_ok"] is False
    assert result["target_anchor_ok"] is False
    assert result["overall_pass"] is False


def test_evaluate_pass_criteria_requires_all_three_conditions():
    """三个条件里只要有一个不达标，overall_pass就应该是False——不是"多数通过就行"。"""
    only_rmsd_bad = MDStabilityResult(
        "DEMO-PARTIAL", ligand_rmsd_last_15ns_mean_angstrom=4.1, hinge_hbond_occupancy_pct=85,
        target_anchor_occupancy_pct=55, md_length_ns=20, n_replicas=3,
    )
    result = evaluate_pass_criteria(only_rmsd_bad)
    assert result["rmsd_ok"] is False
    assert result["hinge_hbond_ok"] is True
    assert result["target_anchor_ok"] is True
    assert result["overall_pass"] is False


def test_run_protein_ligand_complex_md_honestly_fails_without_openmm_stack():
    """
    主.venv(Python 3.9)没有装openmm/openff-toolkit/openmmforcefields(真实调用
    需要一个独立的conda环境，2026-10起已经真实装好并验证跑通，见
    doc/04-funnel-docking.md)。同样的"不是假设没装"原则：这个测试检查的是
    lazy import失败时的诚实返回。
    """
    result = run_protein_ligand_complex_md(
        ligand_id="DEMO-001",
        receptor_pdb_path="/nonexistent/receptor.pdb",
        ligand_sdf_path="/nonexistent/ligand.sdf",
        out_trajectory_path="/tmp/should_not_be_created.dcd",
    )
    if not result["ok"]:
        assert "openmm" in result["reason"] or "openff" in result["reason"]
    else:
        raise AssertionError("不存在的文件路径不应该成功返回ok=True")
