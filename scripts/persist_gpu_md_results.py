"""
把`notebooks/colab_gpu_md.ipynb`第10步真实跑出来的`my_candidates_md_results.json`
(真实候选分子在免费GPU上跑出来的复合物MD结果)落库，接上第1步
`scripts/export_top_candidates_for_gpu_md.py`导出的那批候选——这是环节4
L3(MD稳定性)第一次有真实候选分子(不是demo分子)的结果进数据库。

**老实说清楚这次落库和`core/md_stability.py::evaluate_pass_criteria()`的差距**：
那个函数要求同时有rmsd/铰链氢键占有率/锚点残基占有率三项才能判定"是否通过"，
但`run_protein_ligand_complex_md()`(本地和notebook里用的是同一套逻辑)目前只
真实算出配体RMSD，氢键占有率两项还没实现(需要按残基名追踪距离，是已知、
写在模块docstring里的差距)。这个脚本不会假装调用`evaluate_pass_criteria()`
拿到一个"通过/不通过"的判定——只老实落库真实算出来的RMSD均值，`qc_pass`
留`NULL`，不编造另外两项数字、也不编造一个基于不完整信息的通过/不通过结论。

**2026-10新增：自动给"持续漂移"的轨迹打假阳性归因标签**——用户要求"后续是不是
要记录这种形态不行，下次生成药物分子就避免"。这里用一个很朴素的启发式(前半段
均值 vs 后半段均值，差超过0.5Å才算有明显趋势)把轨迹分成settling/drifting/
oscillating三类；只有"drifting"(后半段持续比前半段差，姿态在往外飘)这一类
会真实写一条`false_positive_attributions`记录，标签是`pose_error`(这7个
枚举值里语义最贴切的一个：L1对接给出的姿态，被后续MD证据认为不是真实稳定的
结合模式)，`attributed_by`老实标成这个启发式的名字，不冒充人工复盘。

**特别要老实说明的边界**：`core/feedback.py::rule_precision_tracking()`的
治理原则要求precision>0.8且≥3个独立化合物验证过，才能把一条规则升级成硬
拒绝——这个脚本只是在**记录**假阳性归因，远没有到"可以变成一条新的
structural_alerts规则去拒绝某个R基团"这一步。目前只有1个真实的drifting
样本(RTA-0009)，不够格做任何统计意义上的归因，这次落库的价值是**开始积累
数据**，不是宣称已经找到了原因——详见`doc/routes/route-a-shp2-sos1.md`里
对R1基团的假设性讨论(明确标注是假设，不是结论)。

用法：
    .venv/bin/python scripts/persist_gpu_md_results.py \
        notebooks/my_candidates_gpu_md_results/my_candidates_md_results.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from db.connect import get_connection
from db.repository import upsert_false_positive_attribution, upsert_funnel_score, upsert_md_stability_qc


def classify_trend(trace: list[float], threshold_angstrom: float = 0.5) -> str:
    """前半段均值 vs 后半段均值，差超过阈值才判定有明显趋势，否则算震荡。"""
    half = len(trace) // 2
    first_half_mean = sum(trace[:half]) / half
    second_half_mean = sum(trace[half:]) / (len(trace) - half)
    if first_half_mean - second_half_mean > threshold_angstrom:
        return "settling"
    if second_half_mean - first_half_mean > threshold_angstrom:
        return "drifting"
    return "oscillating"


def persist(results_path: Path) -> None:
    with open(results_path) as f:
        results = json.load(f)

    conn = get_connection()
    for compound_id, result in results.items():
        if not result.get("ok"):
            print(f"[跳过] {compound_id}: 这次GPU MD本身没有成功({result.get('error', '原因未知')})")
            continue

        cur = conn.cursor()
        cur.execute(
            "SELECT id, ensemble_id FROM funnel_scores WHERE compound_id=%s AND funnel_level='L1' LIMIT 1",
            (compound_id,),
        )
        row = cur.fetchone()
        if row is None:
            print(f"[跳过] {compound_id}: 数据库里没有这个候选的L1记录，不知道该关联哪个ensemble_id")
            continue
        l1_funnel_score_id, ensemble_id = row

        trace = result["ligand_rmsd_trace_angstrom"]
        mean_rmsd = round(sum(trace) / len(trace), 2)

        upsert_funnel_score(
            conn,
            compound_id=compound_id,
            ensemble_id=ensemble_id,
            funnel_level="L3",
            method="md_stability",
            score_value=None,  # 没有算MM-GBSA能量，这一级目前只有稳定性指标，不编造一个kcal/mol数字
            n_poses_or_frames=len(trace),
            passed_filter=None,  # 三项pass criteria缺两项，不能下通过/不通过的结论，老实留空
            raw_output_path=f"notebooks/my_candidates_gpu_md_results/{compound_id}_md.dcd",
        )
        cur.execute(
            "SELECT id FROM funnel_scores WHERE compound_id=%s AND ensemble_id=%s AND funnel_level='L3'",
            (compound_id, ensemble_id),
        )
        funnel_score_id = cur.fetchone()[0]

        upsert_md_stability_qc(
            conn,
            funnel_score_id=funnel_score_id,
            rmsd_angstrom=mean_rmsd,
            hinge_hbond_occupancy_pct=None,
            target_anchor_occupancy_pct=None,
            qc_pass=None,
        )
        trend = classify_trend(trace)
        print(f"[落库] {compound_id}: 真实GPU({result.get('platform')})跑出均值RMSD={mean_rmsd}Å，"
              f"轨迹趋势={trend}，已写入funnel_scores(id={funnel_score_id})+md_stability_qc")

        if trend == "drifting":
            half = len(trace) // 2
            first_mean = round(sum(trace[:half]) / half, 2)
            second_mean = round(sum(trace[half:]) / (len(trace) - half), 2)
            upsert_false_positive_attribution(
                conn,
                compound_id=compound_id,
                funnel_score_id=l1_funnel_score_id,
                attribution_label="pose_error",
                attributed_by="automated_rmsd_trend_heuristic_v1",
                notes=(
                    f"L3 MD真实轨迹持续漂移(前半段均值{first_mean}Å -> 后半段均值{second_mean}Å)，"
                    f"L1对接给出的姿态可能不是真实稳定的结合模式。轨迹只有{len(trace)}帧/"
                    f"{result.get('n_steps_run', '?')}步(远短于真实15ns窗口要求)，这是单个候选的"
                    "自动化初筛信号，不是人工复盘结论，也还不构成能升级成structural_alerts规则的"
                    "统计证据(需要precision>0.8且≥3个独立化合物，见core/feedback.py治理原则)。"
                ),
            )
            print(f"       -> 轨迹持续漂移，已记录false_positive_attribution(pose_error)，"
                  f"供后续积累更多样本后做规则治理分析")

    conn.commit()


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(1)
    persist(Path(sys.argv[1]))
