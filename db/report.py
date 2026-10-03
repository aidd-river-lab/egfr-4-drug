"""
把数据库里一轮计算结果整理成人读的报告。查数字只是第一步，这个模块负责把
"这些数字意味着什么"讲清楚——术语对照见 doc/learning-guide.md。

用法：
    .venv/bin/python -m db.report                      # 三条路线都打印
    .venv/bin/python -m db.report --route route_c_4th_gen_tki
"""
from __future__ import annotations

import argparse

from db.connect import get_connection

ROUTE_LABELS = {
    "route_a_shp2_sos1": "路线A：SHP2/SOS1变构抑制剂",
    "route_b_degrader": "路线B：突变选择性EGFR降解剂",
    "route_c_4th_gen_tki": "路线C：四代EGFR TKI(C797S)",
}


def _fetch_dicts(cur, sql, params=()):
    cur.execute(sql, params)
    cols = [d[0] for d in cur.description]
    return [dict(zip(cols, row)) for row in cur.fetchall()]


def print_route_report(conn, route_id: str) -> None:
    cur = conn.cursor()
    label = ROUTE_LABELS.get(route_id, route_id)
    print(f"\n{'=' * 70}\n{label}\n{'=' * 70}")

    cur.execute("SELECT COUNT(*) FROM compounds WHERE route_id = %s", (route_id,))
    n_compounds = cur.fetchone()[0]
    if n_compounds == 0:
        print("数据库里还没有这条路线的数据——先跑对应workflow并传persist=True。")
        return
    print(f"候选分子总数: {n_compounds}")

    cur.execute(
        "SELECT COUNT(*) FROM structural_alert_hits h JOIN compounds c ON h.compound_id = c.compound_id "
        "WHERE c.route_id = %s AND h.action = 'REJECT'",
        (route_id,),
    )
    n_rejected = cur.fetchone()[0]
    print(f"结构警示REJECT淘汰: {n_rejected} 个"
          + ("（这批全部通过，说明R基团库本身选得比较保守/干净）" if n_rejected == 0 else ""))

    # ADMET画像统计
    rows = _fetch_dicts(
        cur,
        "SELECT a.* FROM admet_descriptors a JOIN compounds c ON a.compound_id = c.compound_id WHERE c.route_id = %s",
        (route_id,),
    )
    if rows:
        mws = [r["mw"] for r in rows if r["mw"] is not None]
        cns = [r["cns_mpo"] for r in rows if r["cns_mpo"] is not None]
        print(f"\nADMET画像: MW范围 {min(mws):.0f}-{max(mws):.0f} Da，CNS MPO范围 {min(cns):.2f}-{max(cns):.2f}"
              f"（阈值4.0以上算及格；本仓库的CNS MPO是近似值，精度边界见doc/06-admet.md）")
        n_pass_cns = sum(1 for c in cns if c >= 4.0)
        print(f"CNS MPO >= 4.0 的候选: {n_pass_cns}/{len(cns)}")

    # 真实L1对接分数(2026-10起可用，见scripts/setup_docking_env.sh)
    docked = _fetch_dicts(
        cur,
        "SELECT f.compound_id, f.ensemble_id, f.score_value FROM funnel_scores f "
        "JOIN compounds c ON f.compound_id = c.compound_id WHERE c.route_id = %s AND f.funnel_level = 'L1' "
        "ORDER BY f.score_value",
        (route_id,),
    )
    if docked:
        scores = [d["score_value"] for d in docked]
        print(f"\n真实L1对接(AutoDock Vina, ensemble={docked[0]['ensemble_id']}): "
              f"{len(docked)}个候选，结合能范围 {min(scores):.2f} ~ {max(scores):.2f} kcal/mol")
        print(f"  最佳: {docked[0]['compound_id']} ({docked[0]['score_value']:.2f})， "
              f"最弱: {docked[-1]['compound_id']} ({docked[-1]['score_value']:.2f})")
    else:
        print("\n还没有真实L1对接数据——用 scripts/setup_docking_env.sh 搭环境后跑 scripts/batch_dock.py")

    # 决策批次(路线A/C才有；路线B目前没有决策阶段)
    rounds = _fetch_dicts(
        cur, "SELECT * FROM decision_rounds WHERE round_id LIKE %s ORDER BY round_date DESC LIMIT 1", (f"{route_id}%",)
    )
    if not rounds:
        print("\n这条路线还没有跑过决策阶段(环节8)——见该路线workflow文件docstring的说明。")
        cur.close()
        return

    round_id = rounds[0]["round_id"]
    print(f"\n最近一轮决策(round_id={round_id}, batch_size={rounds[0]['batch_size']}):")

    members = _fetch_dicts(
        cur,
        """
        SELECT m.compound_id, m.batch_role, m.is_pareto_optimal, m.desirability_score,
               a.mw, a.cns_mpo, s.sa_score, f.score_value AS dock_score
        FROM decision_batch_members m
        JOIN admet_descriptors a ON m.compound_id = a.compound_id
        LEFT JOIN synthesis_records s ON m.compound_id = s.compound_id
        LEFT JOIN funnel_scores f ON m.compound_id = f.compound_id AND f.funnel_level = 'L1'
        WHERE m.round_id = %s
        ORDER BY FIELD(m.batch_role, 'exploit', 'explore', 'hypothesis_test', 'control'), m.desirability_score DESC
        """,
        (round_id,),
    )
    role_explain = {
        "exploit": "exploit(利用)：当前模型认为最好的——desirability分数(CNS MPO和SA score的几何平均)最高",
        "explore": "explore(探索)：按不确定度选的，本轮没有不确定度数据，诚实留空",
        "hypothesis_test": "hypothesis_test(假设检验)：药化指定的对照实验，本轮没有人工指定，诚实留空",
        "control": "control(对照)：从被结构警示/低分淘汰的候选池里随机抽的，检验规则有没有误杀好分子",
    }
    seen_roles = set()
    for m in members:
        role = m["batch_role"]
        if role not in seen_roles:
            print(f"\n  [{role_explain.get(role, role)}]")
            seen_roles.add(role)
        pareto_tag = "Pareto最优" if m["is_pareto_optimal"] else ""
        dock_str = f"真实对接={m['dock_score']:.2f}" if m["dock_score"] is not None else "无对接数据"
        print(
            f"    {m['compound_id']}: desirability={m['desirability_score']:.3f} "
            f"MW={m['mw']:.0f} CNS_MPO={m['cns_mpo']:.2f} SA_score={m['sa_score']:.2f} {dock_str} {pareto_tag}"
        )
    cur.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--route", choices=list(ROUTE_LABELS.keys()), default=None)
    args = parser.parse_args()

    conn = get_connection()
    try:
        routes = [args.route] if args.route else list(ROUTE_LABELS.keys())
        for route_id in routes:
            print_route_report(conn, route_id)
    finally:
        conn.close()
