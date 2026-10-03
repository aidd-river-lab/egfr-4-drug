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

from core.route_config import route_config_dir


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


def evaluate_pass_criteria(result: MDStabilityResult, config_dir: Path | None = None) -> dict:
    """
    按 pipeline.yaml 的 funnel.L3.pass_criteria 判定这个MD稳定性结果是否通过。
    这是真实、可测试的判定逻辑，哪怕 result 本身来自占位/演示数据也能验证逻辑对不对。

    config_dir: 不传则默认路线C的pipeline.yaml；路线A/B传各自的 route_config_dir(...)。
    注意"target_anchor"这个key名字历史上是照着路线C的Ser797场景起的，路线A/B复用这个
    判定函数时，可以在自己的pipeline.yaml里把这个key挪用来表示别的锚点概念(比如路线A
    SHP2 tunnel位点的关键残基接触)，不需要改这段代码。
    """
    config_dir = config_dir or route_config_dir()
    with open(config_dir / "pipeline.yaml", encoding="utf-8") as f:
        criteria = yaml.safe_load(f)["funnel"]["L3"]["pass_criteria"]

    rmsd_ok = result.ligand_rmsd_last_15ns_mean_angstrom < criteria["ligand_rmsd_last_15ns_mean_max_angstrom"]
    hinge_ok = result.hinge_hbond_occupancy_pct > criteria["hinge_hbond_occupancy_min_pct"]
    anchor_ok = result.target_anchor_occupancy_pct > criteria["target_anchor_occupancy_min_pct"]
    # anchor具体指什么(Ser797-OG还是别的路线的锚点残基)是路线专属知识，不写死在代码里，
    # 从各自pipeline.yaml的target_anchor_note读，没配就给个通用占位说明
    anchor_note = criteria.get("target_anchor_note", "target_anchor锚点残基含义见该路线pipeline.yaml")

    return {
        "ligand_id": result.ligand_id,
        "rmsd_ok": rmsd_ok,
        "hinge_hbond_ok": hinge_ok,
        "target_anchor_ok": anchor_ok,
        "overall_pass": rmsd_ok and hinge_ok and anchor_ok,
        "note": f"{anchor_note}；这里仍按三条件都满足才算L3通过，是否要把anchor降级成加分项"
        "由调用方(core/decide.py)在更上层的多目标打分里处理，不在这一级漏判",
    }


if __name__ == "__main__":
    print("=== run_md_stability_stub (预期 ok=False) ===")
    print(run_md_stability_stub("DEMO-001", "receptor_clean.pdb"))

    print("\n=== evaluate_pass_criteria (用合成数据验证判定逻辑本身是对的) ===")
    good = MDStabilityResult("DEMO-GOOD", ligand_rmsd_last_15ns_mean_angstrom=1.2, hinge_hbond_occupancy_pct=85, target_anchor_occupancy_pct=55, md_length_ns=20, n_replicas=3)
    bad = MDStabilityResult("DEMO-BAD", ligand_rmsd_last_15ns_mean_angstrom=4.1, hinge_hbond_occupancy_pct=30, target_anchor_occupancy_pct=10, md_length_ns=20, n_replicas=3)
    print(evaluate_pass_criteria(good))
    print(evaluate_pass_criteria(bad))
