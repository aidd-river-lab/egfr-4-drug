"""
端到端工作流：路线B(双功能降解剂)的 warhead x linker x E3配体 组合枚举 + ADMET画像。

!! 范围声明，比 run_discovery_round.py 更有限 !!
这个工作流能跑的部分：环节3(三组分枚举)/环节1(标准化)/环节6(ADMET描述符+结构警示)/
环节7.1(SA score)。

它**不**包含多目标决策(环节8)的Pareto/批次选择——因为route_b的核心TPP指标
(dc50_nm/dmax_pct/hook_effect_margin)都需要warhead/E3配体各自的结合常数(Kd)
作为输入，而这些Kd要么来自真实FEP/对接(环节4，本仓库诚实占位)，要么来自已发表的
E3配体文献值(可以查到，但本工作流没有接入)。**没有真实Kd就不编一个假的去跑决策
算法**——这和项目贯穿始终的原则一致(宁可少做一步，也不让人误以为Pareto前沿是
基于真实数据算出来的)。core/ternary_complex.py::evaluate_hook_effect_risk()的
演示见该模块自己的__main__。

用法：
    .venv/bin/python -m workflows.run_bifunctional_round                 # 不落库
    .venv/bin/python -c "from workflows.run_bifunctional_round import run; run(persist=True)"  # 需要.env
"""
from __future__ import annotations

import pandas as pd

from core.admet import compute_cns_mpo, compute_descriptors, estimate_herg_risk_proxy, run_structural_alerts
from core.enumerate import enumerate_bifunctional_library
from core.route_config import route_config_dir
from core.standardize import standardize_batch
from core.synthesis import compute_sa_score

ROUTE_ID = "route_b_degrader"


def run(max_combinations: int | None = None, persist: bool = False) -> pd.DataFrame:
    config_dir = route_config_dir(ROUTE_ID)
    report_lines = [f"[路线] {ROUTE_ID}"]

    # ---- 环节3: warhead x linker x E3配体 三组分枚举 ----
    enumerated = enumerate_bifunctional_library(config_dir=config_dir, max_combinations=max_combinations)
    report_lines.append(f"[环节3 三组分枚举] -> {len(enumerated)} 个连通、去重后的双功能分子候选")
    r_group_cols = [c for c in enumerated.columns if c.endswith("_name")]

    # ---- 环节1: 标准化 + 去重(InChIKey) ----
    std_result = standardize_batch(enumerated["smiles"].tolist())
    report_lines.append(
        f"[环节1 标准化] {len(std_result.curated)} 成功 / {len(std_result.quarantine)} 隔离"
        f"(success_rate={std_result.success_rate:.1%})"
    )
    curated = std_result.curated.drop_duplicates(subset="inchikey").reset_index(drop=True)
    curated["compound_id"] = [f"RTB-{i:04d}" for i in range(len(curated))]
    curated = curated.merge(enumerated[["smiles", *r_group_cols]], left_on="raw_smiles", right_on="smiles", how="left")

    # ---- 环节6: ADMET描述符 + 结构警示 ----
    # 注意：compute_cns_mpo()的desirability曲线是按传统小分子(MW<500)校准的，双功能分子
    # 天然beyond Rule-of-5(本次demo MW普遍700+)，算出来的cns_mpo分数会系统性偏低，
    # 这是预期的，不代表这些分子"不行"——只是说明这把尺子不是为这类分子设计的，
    # 真实决策需要专门针对beyond-Ro5分子校准的CNS渗透模型，这里先如实展示数字。
    rows = []
    alert_hit_rows = []
    for _, row in curated.iterrows():
        smiles = row["standardized_smiles"]
        desc = compute_descriptors(smiles)
        mpo = compute_cns_mpo(desc) if desc.get("ok") else {"ok": False, "cns_mpo": None, "pass_threshold_4_0": None}
        herg = estimate_herg_risk_proxy(smiles)
        alerts = run_structural_alerts(smiles)
        sa = compute_sa_score(smiles)
        for h in alerts.hits:
            alert_hit_rows.append((row["compound_id"], h.rule_name, h.severity, h.action, h.advice))
        rows.append(
            {
                "compound_id": row["compound_id"],
                "smiles": smiles,
                "mw": desc.get("mw"),
                "clogp": desc.get("clogp"),
                "clogd_approx": desc.get("clogd_approx"),
                "tpsa": desc.get("tpsa"),
                "hbd": desc.get("hbd"),
                "hba": desc.get("hba"),
                "fsp3": desc.get("fsp3"),
                "has_basic_aliphatic_amine": desc.get("has_basic_aliphatic_amine"),
                "pka_proxy": desc.get("pka_proxy"),
                "cns_mpo": mpo.get("cns_mpo"),  # 传统小分子尺度，见本函数开头注释
                "cns_mpo_pass_4_0": mpo.get("pass_threshold_4_0"),
                "herg_risk_proxy": herg.get("herg_risk_proxy"),
                "sa_score": sa.get("sa_score"),
                "alert_verdict": alerts.verdict,
            }
        )
    candidates = pd.DataFrame(rows)
    report_lines.append(f"[环节6 ADMET画像] 完成，{len(candidates)} 个候选(注意cns_mpo是传统小分子尺度，见代码注释)")

    print("\n".join(report_lines))
    print()
    print("=== 候选分子ADMET画像(没有Pareto/批次选择——原因见本文件docstring) ===")
    print(candidates.to_string())

    if persist:
        from db.connect import get_connection
        from db.repository import persist_discovery_round

        conn = get_connection()
        try:
            summary = persist_discovery_round(
                conn,
                route_id=ROUTE_ID,
                # designed_molecules.scaffold_id是NOT NULL——路线B没有真正的"骨架"概念，
                # 这里用一个固定的分类标签占位，具体的warhead/linker/e3_ligand组合在
                # r_group_assignment JSON里，不在这个字段里
                scaffold_id="bifunctional",
                enumeration_batch=f"bifunctional-{len(enumerated)}",
                r_group_cols=r_group_cols,
                curated_df=curated,
                quarantine_df=std_result.quarantine,
                all_candidates_df=candidates,
                alert_hit_rows=alert_hit_rows,
                batch_combined_df=None,  # 路线B还没有决策阶段，见本文件docstring
            )
            print("\n=== 已落库(MySQL) ===")
            print(summary)
        finally:
            conn.close()

    return candidates


if __name__ == "__main__":
    run()

# TODO(有真实Kd数据或三元复合物预测工具后按这里接入):
#   - 从文献查lenalidomide/pomalidomide对CRBN的真实Kd，作为e3_ligand_kd_nm输入
#   - warhead_kd_nm 需要环节4 L1-L2对接或FEP产出(本仓库诚实占位，接入后才有真实值)
#   - 有了两个Kd后，用 core/ternary_complex.py::evaluate_hook_effect_risk() 评估每个候选
#   - 三元复合物cooperativity_alpha/DC50/Dmax需要PRosettaC/Rosetta类工具，本仓库未安装
#   - 真实数据齐全后，才能用 core/decide.py 对 dc50_nm/dmax_pct/hook_effect_margin 跑Pareto/批次选择
