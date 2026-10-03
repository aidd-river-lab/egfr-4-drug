"""
端到端工作流：跑一轮"纯RDKit可计算"的DMTA循环。路线A/C共用这个工作流
(都是固定骨架+R基团的枚举拓扑)，路线B用 run_bifunctional_round.py。

!! 范围声明 !!
这个工作流串联的是环节1(标准化)/3(枚举)/6(ADMET规则)/7.1(SA score)/8(决策)——
也就是不需要GPU、MD引擎、训练好的ML模型就能在任何机器上真实跑通的部分。

它**不**包含：
  - 环节2 结构系综准备(需要PDB/Rosetta/OpenMM)
  - 环节4 L1-L4 对接/MM-GBSA/FEP打分漏斗(需要Vina+Meeko或GPU集群)
  - 环节4.3 共价对接(需要CovDock/AutoDock4-covalent)
  - 环节5 选择性引擎(依赖环节4的FEP map)
  - 环节9 湿实验回传分析(需要真实湿实验数据)
这些环节的真实计算函数都已经在对应的 core/*.py 里实现为正确接口+诚实占位
(调用会返回 ok=False 和具体缺什么基础设施，不会编造数值)，可以在有算力/数据后接入
这个工作流，位置见本文件最后的 TODO 注释。

用法：
    .venv/bin/python -m workflows.run_discovery_round                          # 路线C(默认)，不落库
    .venv/bin/python -c "from workflows.run_discovery_round import run; run(route_id='route_a_shp2_sos1', scaffold_id='demo_aminopyrazine_tunnel')"
    .venv/bin/python -c "from workflows.run_discovery_round import run; run(persist=True)"  # 需要项目根目录下有.env(见.env.example)
"""
from __future__ import annotations

import pandas as pd

from core.admet import compute_cns_mpo, compute_descriptors, estimate_herg_risk_proxy, run_structural_alerts
from core.decide import DesirabilitySpec, desirability_score, pareto_front, select_batch
from core.enumerate import enumerate_from_scaffold
from core.route_config import route_config_dir
from core.standardize import standardize_batch
from core.synthesis import compute_sa_score

DEFAULT_SCAFFOLD_BY_ROUTE = {
    "route_a_shp2_sos1": "demo_aminopyrazine_tunnel",
    "route_c_4th_gen_tki": "demo_aminopyrimidine_biphenyl",
}

# compound_id前缀按route_id走，不能按scaffold_id截断生成——这里曾经真实踩过一个坑：
# "demo_aminopyrazine_tunnel"[:12] 和 "demo_aminopyrimidine_biphenyl"[:12] 都是
# "demo_aminopy"，导致路线A和路线C生成了同一批compound_id("demo_aminopy-0000"等)，
# 落库时两条路线的分子在compounds表里(PRIMARY KEY只有compound_id)互相覆盖，
# route_id字段显示的是先插入那条路线的、但smiles_canonical却是后插入那条路线的——
# 一次真实的跨路线数据污染事故。现在按route_id生成前缀，不同路线永远不会撞ID。
COMPOUND_ID_PREFIX_BY_ROUTE = {
    "route_a_shp2_sos1": "RTA",
    "route_c_4th_gen_tki": "RTC",
}


def run(
    route_id: str = "route_c_4th_gen_tki",
    scaffold_id: str | None = None,
    batch_size: int = 20,
    persist: bool = False,
) -> pd.DataFrame:
    scaffold_id = scaffold_id or DEFAULT_SCAFFOLD_BY_ROUTE[route_id]
    config_dir = route_config_dir(route_id)
    report_lines = [f"[路线] {route_id}"]

    # ---- 环节3: 枚举 ----
    enumerated = enumerate_from_scaffold(scaffold_id, config_dir=config_dir)
    report_lines.append(f"[环节3 枚举] {scaffold_id} -> {len(enumerated)} 个连通、去重后的候选SMILES")
    r_group_cols = [c for c in enumerated.columns if c.endswith("_name")]

    # ---- 环节1: 标准化 + 去重(InChIKey) ----
    std_result = standardize_batch(enumerated["smiles"].tolist())
    report_lines.append(
        f"[环节1 标准化] {len(std_result.curated)} 成功 / {len(std_result.quarantine)} 隔离"
        f"(success_rate={std_result.success_rate:.1%})"
    )
    curated = std_result.curated.drop_duplicates(subset="inchikey").reset_index(drop=True)
    id_prefix = COMPOUND_ID_PREFIX_BY_ROUTE[route_id]
    curated["compound_id"] = [f"{id_prefix}-{i:04d}" for i in range(len(curated))]
    # 带上R基团信息，供 db/repository.py::upsert_designed_molecules() 落库用
    curated = curated.merge(enumerated[["smiles", *r_group_cols]], left_on="raw_smiles", right_on="smiles", how="left")

    # ---- 环节6: ADMET描述符 + CNS MPO + 结构警示 ----
    # 注意：这一步对curated里**全部**候选计算(包括后面会被REJECT门槛淘汰的)，
    # 因为落库时这些数据本身是有价值的记录，不应该因为没进最终批次就被丢弃
    # (见 db/repository.py::persist_discovery_round 的docstring)。
    rows = []
    alert_hit_rows = []  # (compound_id, rule_name, severity, action, advice)，落库给 structural_alert_hits
    for _, row in curated.iterrows():
        smiles = row["standardized_smiles"]
        desc = compute_descriptors(smiles)
        mpo = compute_cns_mpo(desc) if desc.get("ok") else {"ok": False, "cns_mpo": None, "pass_threshold_4_0": None}
        herg = estimate_herg_risk_proxy(smiles)
        alerts = run_structural_alerts(smiles)  # 这份demo骨架不走共价路线，不传warhead_smarts
        sa = compute_sa_score(smiles)
        for h in alerts.hits:
            alert_hit_rows.append((row["compound_id"], h.rule_name, h.severity, h.action, h.advice))
        rows.append(
            {
                "compound_id": row["compound_id"],
                "smiles": smiles,
                "inchikey": row["inchikey"],
                "mw": desc.get("mw"),
                "clogp": desc.get("clogp"),
                "clogd_approx": desc.get("clogd_approx"),
                "tpsa": desc.get("tpsa"),
                "hbd": desc.get("hbd"),
                "hba": desc.get("hba"),
                "fsp3": desc.get("fsp3"),
                "has_basic_aliphatic_amine": desc.get("has_basic_aliphatic_amine"),
                "pka_proxy": desc.get("pka_proxy"),
                "cns_mpo": mpo.get("cns_mpo"),
                "cns_mpo_pass_4_0": mpo.get("pass_threshold_4_0"),
                "herg_risk_proxy": herg.get("herg_risk_proxy"),
                "sa_score": sa.get("sa_score"),
                "alert_verdict": alerts.verdict,
                "alert_hits": ",".join(h.rule_name for h in alerts.hits) or "none",
            }
        )
    all_candidates = pd.DataFrame(rows)

    n_before_gate = len(all_candidates)
    rejected = all_candidates[all_candidates["alert_verdict"] == "REJECT"].copy()
    candidates = all_candidates[all_candidates["alert_verdict"] != "REJECT"].reset_index(drop=True)
    report_lines.append(
        f"[环节6 结构警示硬门] {n_before_gate} -> {len(candidates)} 通过(REJECT淘汰 {len(rejected)} 个)"
    )

    # ---- 环节8: 多目标决策 ----
    pf = pareto_front(candidates, {"cns_mpo": "higher_is_better", "sa_score": "lower_is_better"})
    candidates["is_pareto_optimal"] = pf
    report_lines.append(f"[环节8 Pareto前沿] {int(pf.sum())} / {len(candidates)} 个候选在前沿上")

    specs = [
        DesirabilitySpec("cns_mpo", lambda x: 1.0 if x >= 4.0 else max(0.0, x / 4.0), weight=1.0),
        DesirabilitySpec("sa_score", lambda x: 1.0 if x <= 3.0 else max(0.0, 1 - (x - 3.0) / 3.0), weight=1.0),
    ]
    candidates["desirability"] = desirability_score(candidates, specs)

    # 诚实声明：这里没有不确定度列(真实不确定度要来自ML ensemble方差或FEP误差，
    # 本轮没有跑那些)，也没有药化指定的假设检验index，所以explore/hypothesis_test
    # 这两档会诚实地留空，不假装有数据去填充(见core/decide.py select_batch的docstring)。
    batch = select_batch(
        candidates, desirability_col="desirability", batch_size=batch_size, rejected_pool=rejected
    )
    report_lines.append(
        f"[环节8 批次选择, batch_size={batch_size}] "
        f"exploit={len(batch.exploit)} explore={len(batch.explore)}(无不确定度数据，诚实留空) "
        f"hypothesis_test={len(batch.hypothesis_test)}(无药化指定假设，诚实留空) "
        f"control={len(batch.control)}(从REJECT池抽样，用于假阴性复活)"
    )

    print("\n".join(report_lines))
    print()
    print("=== 最终批次(combined) ===")
    print(batch.combined[["compound_id", "batch_role", "cns_mpo", "sa_score", "desirability"]].to_string())

    if persist:
        from db.connect import get_connection
        from db.repository import persist_discovery_round

        conn = get_connection()
        try:
            summary = persist_discovery_round(
                conn,
                route_id=route_id,
                scaffold_id=scaffold_id,
                enumeration_batch=f"{scaffold_id}-{len(enumerated)}",
                r_group_cols=r_group_cols,
                curated_df=curated,
                quarantine_df=std_result.quarantine,
                all_candidates_df=all_candidates,
                alert_hit_rows=alert_hit_rows,
                batch_combined_df=batch.combined,
                batch_size=batch_size,
                quota={"exploit_pct": 50, "explore_pct": 25, "hypothesis_test_pct": 15, "control_pct": 10},
            )
            print("\n=== 已落库(MySQL) ===")
            print(summary)
        finally:
            conn.close()

    return batch.combined


if __name__ == "__main__":
    run()

# TODO(有算力/数据后按这里接入，不要改动上面已验证的逻辑):
#   - 环节2: core/structures.py fetch_pdb + mutate_residue_stub -> 产出 StructureEnsembleMember 列表
#   - 环节4 L1-L2: 对每个candidates["smiles"] x 每个ensemble成员跑 core/docking.py run_vina_docking
#   - 环节4 L3-L4: core/md_stability.py / core/mmgbsa.py / core/fep.py，只对L2晋级的子集跑
#   - 环节4.3: 如果chemistry_route是covalent_new_site/covalent_pan_mutant_broad，core/covalent.py evaluate_attack_geometry
#   - 环节5: core/selectivity.py compute_dddg_from_fep_maps，需要环节4.4产出的FEP map
#   - 环节9: 拿到湿实验结果后用 core/feedback.py compare_predicted_vs_actual 复盘本轮
