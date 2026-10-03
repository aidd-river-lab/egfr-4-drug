"""
环节8：多目标决策引擎

设计文档明确反对加权求和(环节8.1)："高分可以靠一项极端值补偿其他项的缺陷"。
这里实现两条推荐路线：
  1. Pareto 前沿非支配排序 —— pareto_front()
  2. 期望效用(geometric mean of desirability) —— desirability_score()，
     几何平均的好处是任何一项接近0会拖垮总分，符合"成药性是乘法不是加法"。

以及环节8.3 的探索/利用/假设检验/对照配额批次选择 —— select_batch()。

本模块是纯算法，不依赖任何化学计算，输入是一个 pandas DataFrame，
可以直接用来演示（见 __main__），也可以接 core/enumerate.py + core/admet.py
算出来的真实候选分子属性表。
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Literal

import numpy as np
import pandas as pd
import yaml

Direction = Literal["higher_is_better", "lower_is_better"]

TARGET_PROFILE_CONFIG = Path(__file__).resolve().parent.parent / "config" / "target_profile.yaml"


def load_decision_config_from_target_profile() -> dict:
    """
    从 target_profile.yaml 读取 decision_engine 配置，对齐 Pareto 目标方向和批次配额。

    target_profile.yaml 里有一条注释："决策引擎(core/decide.py)读取这里的字段名来对齐
    TPP...不要改字段名除非同步改core/decide.py里的TPP_FIELD_MAP"——但此前只有这条注释，
    没有对应代码，配置和代码实际上没有接上。这个函数把这个对接点真正补上。

    decision_engine.pareto_objectives 只给了指标名，方向要从 tpp 表的 direction 字段查。
    tpp 里 direction=in_range 的指标(比如 clogd74)不能直接喂给 pareto_front()(它只认
    单调的 higher/lower_is_better)，这种指标被诚实跳过并列在 skipped_non_monotonic_metrics
    里——调用方需要自己把 in_range 转换成 desirability_score 的一个 spec，而不是指望
    这里偷偷帮你决定怎么转换。
    """
    with open(TARGET_PROFILE_CONFIG, encoding="utf-8") as f:
        profile = yaml.safe_load(f)

    tpp_direction = {row["metric"]: row["direction"] for row in profile["tpp"]}
    requested = profile["decision_engine"]["pareto_objectives"]

    pareto_objectives = {}
    missing_from_tpp = []
    non_monotonic = []
    for metric in requested:
        if metric not in tpp_direction:
            missing_from_tpp.append(metric)
            continue
        direction = tpp_direction[metric]
        if direction not in ("higher_is_better", "lower_is_better"):
            non_monotonic.append(metric)
            continue
        pareto_objectives[metric] = direction

    return {
        "pareto_objectives": pareto_objectives,
        "missing_from_tpp_table": missing_from_tpp,  # 配置写错了/漏填tpp行，属于配置bug，应该被修掉
        "skipped_non_monotonic_metrics": non_monotonic,  # 比如in_range，属于预期行为，不是bug
        "batch_quota": profile["decision_engine"]["batch_quota"],
    }


# ------------------------------------------------------------------
# 1. Pareto 前沿
# ------------------------------------------------------------------
def pareto_front(df: pd.DataFrame, objectives: dict[str, Direction]) -> pd.Series:
    """
    返回一个布尔 Series：True 表示该行在给定目标集合上是非支配解(不存在另一行在所有
    目标上都不比它差、且至少一项严格更好)。

    objectives: {"cell_ic50_nm": "lower_is_better", "wt_selectivity_fold": "higher_is_better", ...}
    """
    if df.empty:
        return pd.Series([], dtype=bool)

    # 统一转成"越大越好"的矩阵，方便比较
    cols = list(objectives.keys())
    values = df[cols].copy()
    for col, direction in objectives.items():
        if direction == "lower_is_better":
            values[col] = -values[col]
    arr = values.to_numpy(dtype=float)

    n = len(arr)
    is_dominated = np.zeros(n, dtype=bool)
    for i in range(n):
        if is_dominated[i]:
            continue
        # j 支配 i 当且仅当 j 在所有目标上 >= i，且至少一项 > i
        ge_all = (arr >= arr[i]).all(axis=1)
        gt_any = (arr > arr[i]).any(axis=1)
        dominators = ge_all & gt_any
        dominators[i] = False
        if dominators.any():
            is_dominated[i] = True

    return pd.Series(~is_dominated, index=df.index, name="is_pareto_optimal")


# ------------------------------------------------------------------
# 2. 期望效用 / MPO 几何平均打分
# ------------------------------------------------------------------
@dataclass
class DesirabilitySpec:
    column: str
    desirability_fn: Callable[[float], float]
    weight: float = 1.0  # 几何平均里的指数权重，默认等权


def desirability_score(df: pd.DataFrame, specs: list[DesirabilitySpec]) -> pd.Series:
    """
    对每个指标算 desirability(0-1)，取加权几何平均作为总分。
    任何一项 desirability=0 会让总分直接变成0（"乘法不是加法"）——这是有意的设计，
    不是bug：一个CNS MPO很好但hERG判为0分的分子，不应该靠其它指标的高分把它捞回来。
    """
    if df.empty:
        return pd.Series([], dtype=float)

    total_weight = sum(s.weight for s in specs)
    log_sum = pd.Series(0.0, index=df.index)
    zero_mask = pd.Series(False, index=df.index)

    for spec in specs:
        d = df[spec.column].apply(spec.desirability_fn).astype(float)
        zero_mask |= d <= 0
        # 用log避免0直接导致log(0)报错，最后统一把zero_mask的位置设成0
        safe_d = d.clip(lower=1e-9)
        log_sum += spec.weight * np.log(safe_d)

    score = np.exp(log_sum / total_weight)
    score[zero_mask] = 0.0
    return score.rename("desirability_score")


# ------------------------------------------------------------------
# 3. 批次选择：利用/探索/假设检验/对照（环节8.3）
# ------------------------------------------------------------------
@dataclass
class BatchSelection:
    exploit: pd.DataFrame
    explore: pd.DataFrame
    hypothesis_test: pd.DataFrame
    control: pd.DataFrame

    @property
    def combined(self) -> pd.DataFrame:
        parts = []
        for role, part in [
            ("exploit", self.exploit),
            ("explore", self.explore),
            ("hypothesis_test", self.hypothesis_test),
            ("control", self.control),
        ]:
            if not part.empty:
                p = part.copy()
                p["batch_role"] = role
                parts.append(p)
        return pd.concat(parts, axis=0) if parts else pd.DataFrame()


def select_batch(
    candidates: pd.DataFrame,
    desirability_col: str,
    batch_size: int,
    uncertainty_col: str | None = None,
    hypothesis_test_indices: list | None = None,
    rejected_pool: pd.DataFrame | None = None,
    quota: dict[str, float] | None = None,
) -> BatchSelection:
    """
    按 pipeline.yaml 里 batch_composition 的配额(默认 50/25/15/10)组一批合成清单。

    hypothesis_test_indices: 这一档本质是药化的化学判断("测试797入口区砜基是否必要"这种
    成对假设)，不是算法能自动生成的——设计文档原话"药化指定"。这里只接受调用方显式传入的
    index 列表，不会凭空生成假设检验分子。如果不传，这一档就是空的，配额会被忽略（不会
    自动挪给其他档位去凑数，免得掩盖"这一轮没做假设检验"这个事实）。

    rejected_pool: 被规则/低分拒绝的候选池，control档从这里随机抽，用来检验规则有没有
    误杀好分子(设计文档环节8.3/9.3反复强调的"假阴性复活机制")。不传就退化为从candidates
    里剩余的随机抽。
    """
    quota = quota or {"exploit_pct": 50, "explore_pct": 25, "hypothesis_test_pct": 15, "control_pct": 10}
    n_exploit = round(batch_size * quota["exploit_pct"] / 100)
    n_explore = round(batch_size * quota["explore_pct"] / 100)
    n_hypothesis = round(batch_size * quota["hypothesis_test_pct"] / 100)
    n_control = batch_size - n_exploit - n_explore - n_hypothesis

    remaining = candidates.copy()

    # exploit: 按desirability取top-N
    exploit = remaining.sort_values(desirability_col, ascending=False).head(n_exploit)
    remaining = remaining.drop(exploit.index)

    # explore: 按不确定度取top-N(没有不确定度列就跳过，不要用假数据填)
    if uncertainty_col and uncertainty_col in remaining.columns:
        explore = remaining.sort_values(uncertainty_col, ascending=False).head(n_explore)
        remaining = remaining.drop(explore.index)
    else:
        explore = remaining.iloc[0:0]

    # hypothesis_test: 只接受显式指定的index
    if hypothesis_test_indices:
        valid_idx = [i for i in hypothesis_test_indices if i in remaining.index]
        hypothesis_test = remaining.loc[valid_idx].head(n_hypothesis)
        remaining = remaining.drop(hypothesis_test.index)
    else:
        hypothesis_test = remaining.iloc[0:0]

    # control: 优先从rejected_pool随机抽（假阴性复活），不够再从剩余candidates里补
    control_parts = []
    if rejected_pool is not None and not rejected_pool.empty:
        n_from_rejected = min(n_control, len(rejected_pool))
        control_parts.append(rejected_pool.sample(n=n_from_rejected, random_state=42))
    n_still_needed = n_control - sum(len(p) for p in control_parts)
    if n_still_needed > 0 and not remaining.empty:
        n_from_remaining = min(n_still_needed, len(remaining))
        control_parts.append(remaining.sample(n=n_from_remaining, random_state=42))
    control = pd.concat(control_parts, axis=0) if control_parts else candidates.iloc[0:0]

    return BatchSelection(exploit=exploit, explore=explore, hypothesis_test=hypothesis_test, control=control)


if __name__ == "__main__":
    print("=== 从 target_profile.yaml 加载决策引擎配置 ===")
    print(load_decision_config_from_target_profile())

    # 合成演示数据(不是真实预测值，只用来演示算法行为)——明确标注避免和真实数据混淆
    rng = np.random.default_rng(42)
    n = 30
    demo = pd.DataFrame(
        {
            "compound_id": [f"DEMO-{i:03d}" for i in range(n)],
            "cell_ic50_nm": rng.uniform(1, 200, n),
            "wt_selectivity_fold": rng.uniform(1, 150, n),
            "cns_mpo": rng.uniform(1, 6, n),
            "uncertainty": rng.uniform(0, 1, n),
        }
    )

    print("=== Pareto front ===")
    pf = pareto_front(
        demo,
        {"cell_ic50_nm": "lower_is_better", "wt_selectivity_fold": "higher_is_better", "cns_mpo": "higher_is_better"},
    )
    print(f"{pf.sum()} / {n} 个在Pareto前沿上")

    print("\n=== Desirability score ===")
    specs = [
        DesirabilitySpec("cell_ic50_nm", lambda x: 1.0 if x <= 20 else max(0.0, 1 - (x - 20) / 80)),
        DesirabilitySpec("wt_selectivity_fold", lambda x: min(1.0, x / 100)),
        DesirabilitySpec("cns_mpo", lambda x: min(1.0, x / 4.5)),
    ]
    demo["desirability"] = desirability_score(demo, specs)
    print(demo.sort_values("desirability", ascending=False).head(5).to_string())

    print("\n=== Batch selection (batch_size=10) ===")
    batch = select_batch(demo, desirability_col="desirability", batch_size=10, uncertainty_col="uncertainty")
    print(f"exploit={len(batch.exploit)} explore={len(batch.explore)} "
          f"hypothesis_test={len(batch.hypothesis_test)} control={len(batch.control)}")
    print(batch.combined[["compound_id", "batch_role", "desirability"]].to_string())
