"""
环节4 L3（前半）：MD 稳定性筛（设计文档称"性价比最高的一步"）

真正需要 OpenMM/GROMACS + GPU 才能产出轨迹，本环境没有GPU，run_md_stability()
诚实返回未执行。但"拿到轨迹统计结果之后怎么判定通过/不通过"这part是纯逻辑，
不需要真的跑过MD就能写对、测对——evaluate_pass_criteria() 是本模块里真正可用、
已经用合成数据验证过的部分。
"""
from __future__ import annotations

import yaml
from dataclasses import dataclass
from pathlib import Path

PIPELINE_CONFIG = Path(__file__).resolve().parent.parent / "config" / "pipeline.yaml"


@dataclass
class MDStabilityResult:
    ligand_id: str
    ligand_rmsd_last_15ns_mean_angstrom: float
    hinge_hbond_occupancy_pct: float
    target_anchor_occupancy_pct: float
    md_length_ns: int
    n_replicas: int
    method: str = "openmm"


def run_md_stability_stub(ligand_id: str, receptor_structure_path: str, md_length_ns: int = 20, n_replicas: int = 3) -> dict:
    """占位：真实运行需要OpenMM + GPU，耗时量级见doc/04-funnel-docking.md的算力预算表。"""
    return {
        "ok": False,
        "reason": f"需要OpenMM+GPU跑{n_replicas}条{md_length_ns}ns轨迹并统计配体RMSD/氢键占据率，本环境无GPU未执行",
        "ligand_id": ligand_id,
    }


def evaluate_pass_criteria(result: MDStabilityResult) -> dict:
    """
    按 pipeline.yaml 的 funnel.L3.pass_criteria 判定这个MD稳定性结果是否通过。
    这是真实、可测试的判定逻辑，哪怕 result 本身来自占位/演示数据也能验证逻辑对不对。
    """
    with open(PIPELINE_CONFIG, encoding="utf-8") as f:
        criteria = yaml.safe_load(f)["funnel"]["L3"]["pass_criteria"]

    rmsd_ok = result.ligand_rmsd_last_15ns_mean_angstrom < criteria["ligand_rmsd_last_15ns_mean_max_angstrom"]
    hinge_ok = result.hinge_hbond_occupancy_pct > criteria["hinge_hbond_occupancy_min_pct"]
    anchor_ok = result.target_anchor_occupancy_pct > criteria["target_anchor_occupancy_min_pct"]

    return {
        "ligand_id": result.ligand_id,
        "rmsd_ok": rmsd_ok,
        "hinge_hbond_ok": hinge_ok,
        "target_anchor_ok": anchor_ok,
        "overall_pass": rmsd_ok and hinge_ok and anchor_ok,
        "note": "target_anchor(Ser797)不是硬门槛，见pipeline.yaml的ser797_hbond_policy，但这里仍按三条件都满足才算L3通过，"
        "是否要把anchor降级成加分项由调用方(core/decide.py)在更上层的多目标打分里处理，不在这一级漏判",
    }


if __name__ == "__main__":
    print("=== run_md_stability_stub (预期 ok=False) ===")
    print(run_md_stability_stub("DEMO-001", "receptor_clean.pdb"))

    print("\n=== evaluate_pass_criteria (用合成数据验证判定逻辑本身是对的) ===")
    good = MDStabilityResult("DEMO-GOOD", ligand_rmsd_last_15ns_mean_angstrom=1.2, hinge_hbond_occupancy_pct=85, target_anchor_occupancy_pct=55, md_length_ns=20, n_replicas=3)
    bad = MDStabilityResult("DEMO-BAD", ligand_rmsd_last_15ns_mean_angstrom=4.1, hinge_hbond_occupancy_pct=30, target_anchor_occupancy_pct=10, md_length_ns=20, n_replicas=3)
    print(evaluate_pass_criteria(good))
    print(evaluate_pass_criteria(bad))
