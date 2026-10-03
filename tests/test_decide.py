"""
core/decide.py 的两条关键设计决策的回归测试：
  1. desirability_score 是几何平均，任何一项为0会把总分拖到0(不是加权和能补偿的)
  2. select_batch 在没有显式传入hypothesis_test_indices时，这一档必须保持空，
     不能为了凑batch_size把配额挪给其他档位（否则会掩盖"这一轮没做假设检验"的事实）
"""
import pandas as pd

from core.decide import (
    DesirabilitySpec,
    desirability_score,
    load_decision_config_from_target_profile,
    pareto_front,
    select_batch,
)


def test_desirability_zero_on_any_dimension_zeros_total_score():
    df = pd.DataFrame({"good_metric": [10.0], "bad_metric": [0.0]})
    specs = [
        DesirabilitySpec("good_metric", lambda x: 1.0),  # 满分
        DesirabilitySpec("bad_metric", lambda x: 0.0),  # 零分
    ]
    score = desirability_score(df, specs)
    assert score.iloc[0] == 0.0


def test_desirability_is_not_a_weighted_sum():
    """一项0.0分和一项1.0分的几何平均必须是0，不能是加权和的0.5。"""
    df = pd.DataFrame({"a": [1.0], "b": [0.0]})
    specs = [DesirabilitySpec("a", lambda x: x, weight=1.0), DesirabilitySpec("b", lambda x: x, weight=1.0)]
    score = desirability_score(df, specs)
    assert score.iloc[0] == 0.0  # 如果是加权和，这里会错误地变成0.5


def test_pareto_front_excludes_strictly_dominated_point():
    df = pd.DataFrame({"x": [1.0, 2.0], "y": [1.0, 2.0]})  # 第0行在两个维度上都被第1行支配
    result = pareto_front(df, {"x": "higher_is_better", "y": "higher_is_better"})
    assert list(result) == [False, True]


def test_select_batch_hypothesis_test_empty_when_no_indices_given():
    candidates = pd.DataFrame({"compound_id": [f"C{i}" for i in range(20)], "desirability": range(20)})
    batch = select_batch(candidates, desirability_col="desirability", batch_size=10)
    assert len(batch.hypothesis_test) == 0


def test_select_batch_does_not_reallocate_empty_hypothesis_quota():
    """
    没传hypothesis_test_indices时，总选出的分子数应该诚实地少于batch_size
    (quota被忽略而不是挪给别的档位去凑数)。
    """
    candidates = pd.DataFrame({"compound_id": [f"C{i}" for i in range(20)], "desirability": range(20)})
    batch = select_batch(candidates, desirability_col="desirability", batch_size=10)
    total = len(batch.exploit) + len(batch.explore) + len(batch.hypothesis_test) + len(batch.control)
    assert total < 10


def test_load_decision_config_resolves_all_pareto_objectives_cleanly():
    """
    回归锁定：decision_engine.pareto_objectives 里列出的每个指标，必须要么在
    pareto_objectives 结果里拿到方向，要么被明确归类成missing_from_tpp_table
    (配置bug)或skipped_non_monotonic_metrics(预期行为，比如in_range指标)。
    目前target_profile.yaml里的四个指标应该全部能在tpp表里找到方向。
    """
    config = load_decision_config_from_target_profile()
    assert config["missing_from_tpp_table"] == []
    assert "cell_ic50_del19_c797s_nm" in config["pareto_objectives"]
    assert "wt_selectivity_fold" in config["pareto_objectives"]
    assert "cns_mpo" in config["pareto_objectives"]
    assert "synthesis_n_steps" in config["pareto_objectives"]


def test_load_decision_config_batch_quota_sums_to_100():
    config = load_decision_config_from_target_profile()
    assert sum(config["batch_quota"].values()) == 100


def test_select_batch_uses_explicit_hypothesis_indices_when_given():
    candidates = pd.DataFrame({"compound_id": [f"C{i}" for i in range(20)], "desirability": range(20)})
    batch = select_batch(
        candidates, desirability_col="desirability", batch_size=10, hypothesis_test_indices=[5, 6, 7]
    )
    assert len(batch.hypothesis_test) > 0
    assert set(batch.hypothesis_test.index).issubset({5, 6, 7})
