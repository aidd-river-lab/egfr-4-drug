"""
环节9：DMTA 闭环、主动学习与 MLOps 的数据层

两部分都是真实可运行的纯逻辑，不需要连真实LIMS就能验证：
  1. WetLabResult：湿实验回传的结构化schema(pydantic强校验，不是自由格式Excel)
  2. compare_predicted_vs_actual / rule_precision_tracking：环节9.2/9.3要求的
     预测-实测复盘和规则精度追踪，纯pandas统计
"""
from __future__ import annotations

from datetime import date
from typing import Literal, Optional

import numpy as np
import pandas as pd
from pydantic import BaseModel, Field
from scipy.stats import spearmanr

FalsePositiveAttribution = Literal[
    "pose_error",  # 位姿错误
    "protonation_state_error",  # 质子化态错误
    "conformer_selection_error",  # 构象选择错误
    "desolvation_underestimated",  # 去溶剂化低估
    "cell_permeability_issue",  # 细胞渗透问题
    "metabolic_instability",  # 代谢不稳定
    "compound_degradation_or_purity",  # 化合物降解或纯度问题
]


class WetLabResult(BaseModel):
    """
    环节9.1 要求的最小字段集，pydantic 在赋值时就强校验类型，不是自由格式Excel。
    smiles_as_made 不能省——合成出来的实际结构和设计稿不一致是家常便饭，
    如果拿实测活性去对设计结构，整个模型会被污染（见环节9.1原文）。
    """

    # 身份
    compound_id: str
    batch_id: str
    smiles_as_made: str
    purity_pct: float = Field(ge=0, le=100)
    chirality_confirmed: bool

    # 活性
    genotype: str
    assay_format: Literal["enzymatic", "cellular", "nanobret"]
    atp_conc_um: Optional[float] = None  # 酶活必填，细胞测定可以是None
    readout: str
    value: float
    unit: str
    censoring: Literal["point", "less_than", "greater_than"] = "point"
    n_replicates: int = Field(ge=1)
    log_sd: Optional[float] = None
    assay_date: date
    plate_id: Optional[str] = None

    # 选择性
    wt_value: Optional[float] = None
    fold_selectivity: Optional[float] = None
    kinase_panel_file: Optional[str] = None

    # ADME
    sol_ph74: Optional[float] = None
    hlm_clint: Optional[float] = None
    rlm_clint: Optional[float] = None
    ppb_pct: Optional[float] = None
    caco2_papp: Optional[float] = None
    mdr1_er: Optional[float] = None
    herg_ic50_um: Optional[float] = None
    gsh_adduct_detected: Optional[bool] = None

    # 体内
    species: Optional[str] = None
    dose_mg_per_kg: Optional[float] = None
    route: Optional[str] = None
    auc: Optional[float] = None
    cmax: Optional[float] = None
    t_half_hours: Optional[float] = None
    f_pct: Optional[float] = None
    kp_uu: Optional[float] = None

    # 元数据
    operator: str
    protocol_version: str
    qc_status: Literal["pass", "fail", "pending"] = "pending"
    comments: str = ""


def compare_predicted_vs_actual(
    df: pd.DataFrame, predicted_col: str, actual_col: str, scaffold_col: str | None = None
) -> dict:
    """
    环节9.2第1条：预测vs实测复盘。整体Spearman rho + MUE，以及(如果提供了scaffold分组)
    按化学型分组的误差——"整体ρ好而某个子系列系统性偏差"是最有价值的信号。
    """
    valid = df.dropna(subset=[predicted_col, actual_col])
    if len(valid) < 3:
        return {"ok": False, "reason": f"有效数据点太少(n={len(valid)})，Spearman相关性至少需要3个点"}

    rho, pvalue = spearmanr(valid[predicted_col], valid[actual_col])
    mue = float(np.mean(np.abs(valid[predicted_col] - valid[actual_col])))

    result = {"ok": True, "n": len(valid), "spearman_rho": round(rho, 3), "p_value": round(pvalue, 4), "mue": round(mue, 3)}

    if scaffold_col and scaffold_col in valid.columns:
        per_scaffold = {}
        for scaffold, group in valid.groupby(scaffold_col):
            if len(group) < 2:
                per_scaffold[scaffold] = {"n": len(group), "mue": None, "note": "样本太少，不单独算MUE"}
                continue
            per_scaffold[scaffold] = {
                "n": len(group),
                "mue": round(float(np.mean(np.abs(group[predicted_col] - group[actual_col]))), 3),
            }
        result["per_scaffold"] = per_scaffold

    return result


def rule_precision_tracking(rule_hits: pd.DataFrame, true_label_col: str, rule_hit_col: str) -> dict:
    """
    环节9.3：规则升级为reject前必须满足 precision > 0.8 且在至少3个独立化合物上验证。
    rule_hits: 每行一个化合物，rule_hit_col=该规则是否命中(bool)，
               true_label_col=事后证实该化合物是否真的有问题(bool，来自湿实验结果)。
    precision = 命中规则且确实有问题 / 命中规则的总数。
    """
    hit = rule_hits[rule_hits[rule_hit_col]]
    if hit.empty:
        return {"ok": False, "reason": "这条规则在提供的数据里一次都没命中，无法算precision"}

    true_positive = hit[true_label_col].sum()
    precision = true_positive / len(hit)
    n_independent_compounds = hit["compound_id"].nunique() if "compound_id" in hit.columns else len(hit)

    eligible_for_reject = precision > 0.8 and n_independent_compounds >= 3
    return {
        "ok": True,
        "n_hits": len(hit),
        "n_true_positive": int(true_positive),
        "precision": round(precision, 3),
        "n_independent_compounds": int(n_independent_compounds),
        "eligible_for_reject_mode": eligible_for_reject,
        "note": "precision>0.8且>=3个独立化合物验证过，才能从warn模式升级成reject(环节9.3硬要求)",
    }


if __name__ == "__main__":
    print("=== WetLabResult schema 校验演示 ===")
    rec = WetLabResult(
        compound_id="INT-0412",
        batch_id="B001",
        smiles_as_made="COc1cc(N(C)CCN(C)C)c(NC(=O)C=C)cc1Nc1nccc(-c2cn(C)c3ccccc23)n1",
        purity_pct=97.5,
        chirality_confirmed=True,
        genotype="del19_C797S",
        assay_format="cellular",
        readout="pIC50",
        value=8.1,
        unit="pIC50",
        n_replicates=3,
        assay_date=date(2026, 10, 1),
        operator="demo",
        protocol_version="v1",
    )
    print(rec.model_dump_json(indent=2)[:400], "...(截断)")

    print("\n=== compare_predicted_vs_actual (合成数据演示) ===")
    rng = np.random.default_rng(0)
    n = 20
    demo_df = pd.DataFrame(
        {
            "predicted_pic50": rng.uniform(6, 9, n),
            "scaffold": rng.choice(["series_A", "series_B"], n),
        }
    )
    demo_df["actual_pic50"] = demo_df["predicted_pic50"] + rng.normal(0, 0.4, n)
    print(compare_predicted_vs_actual(demo_df, "predicted_pic50", "actual_pic50", scaffold_col="scaffold"))

    print("\n=== rule_precision_tracking (合成数据演示) ===")
    rule_df = pd.DataFrame(
        {
            "compound_id": [f"C{i}" for i in range(10)],
            "rule_hit": [True, True, True, True, False, False, True, True, False, True],
            "actually_bad": [True, True, False, True, False, False, True, True, False, True],
        }
    )
    print(rule_precision_tracking(rule_df, true_label_col="actually_bad", rule_hit_col="rule_hit"))
