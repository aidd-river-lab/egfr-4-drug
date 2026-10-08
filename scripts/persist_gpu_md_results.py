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
from db.repository import upsert_funnel_score, upsert_md_stability_qc


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
            "SELECT ensemble_id FROM funnel_scores WHERE compound_id=%s AND funnel_level='L1' LIMIT 1",
            (compound_id,),
        )
        row = cur.fetchone()
        if row is None:
            print(f"[跳过] {compound_id}: 数据库里没有这个候选的L1记录，不知道该关联哪个ensemble_id")
            continue
        ensemble_id = row[0]

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
        trend = "settling(前段高后段低)" if trace[0] > trace[-1] + 0.5 else (
            "drifting(前段低后段高)" if trace[-1] > trace[0] + 0.5 else "oscillating(无明显趋势)"
        )
        print(f"[落库] {compound_id}: 真实GPU({result.get('platform')})跑出均值RMSD={mean_rmsd}Å，"
              f"轨迹趋势={trend}，已写入funnel_scores(id={funnel_score_id})+md_stability_qc")

    conn.commit()


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(1)
    persist(Path(sys.argv[1]))
