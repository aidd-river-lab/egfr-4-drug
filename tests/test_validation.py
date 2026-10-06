"""
core/validation.py 的回归测试。纯数学(AUC的Mann-Whitney U公式、富集因子)，
不需要Vina/PyRosetta，用合成数据验证逻辑本身对不对；真实DUD-E数据跑出来的结果见
validation/dude_egfr/results/ 和 doc/12-retrospective-validation.md。
"""
import random

from core.validation import compute_auc_by_similarity_to_reference, compute_enrichment_factor, compute_roc_auc


def test_auc_perfect_separation_is_one():
    labels = [True, True, True, False, False, False]
    scores = [-9.0, -8.5, -8.0, -5.0, -4.5, -4.0]  # 活性分子分数全部更低(更强结合)
    result = compute_roc_auc(labels, scores)
    assert result["ok"] is True
    assert result["auc"] == 1.0


def test_auc_perfectly_reversed_is_zero():
    labels = [True, True, True, False, False, False]
    scores = [-4.0, -4.5, -5.0, -8.0, -8.5, -9.0]  # 活性分子反而分数更差
    result = compute_roc_auc(labels, scores)
    assert result["auc"] == 0.0


def test_auc_ties_count_as_half():
    """两个正例和负例分数完全相等——按惯例算0.5，不是0也不是1。"""
    labels = [True, False]
    scores = [-7.0, -7.0]
    result = compute_roc_auc(labels, scores)
    assert result["auc"] == 0.5


def test_auc_converges_to_half_for_independent_random_data():
    """标签和分数完全独立生成时，n足够大，AUC应该落在0.5附近(留足够宽的容忍区间，
    这是统计性质不是精确值，不能卡得太死导致测试flaky)。"""
    rng = random.Random(0)
    labels = [rng.random() < 0.5 for _ in range(400)]
    scores = [rng.uniform(-10, -4) for _ in range(400)]
    result = compute_roc_auc(labels, scores)
    assert 0.4 < result["auc"] < 0.6


def test_auc_honestly_fails_without_both_classes():
    result = compute_roc_auc([True, True], [-9.0, -8.0])
    assert result["ok"] is False


def test_auc_honestly_fails_on_length_mismatch():
    result = compute_roc_auc([True, False], [-9.0])
    assert result["ok"] is False


def test_auc_higher_score_is_better_direction_flag():
    """lower_score_is_more_positive=False时，分数越高应该越被判定为活性——
    用和test_auc_perfect_separation_is_one反过来的数据验证方向标志位真的生效。"""
    labels = [True, True, True, False, False, False]
    scores = [9.0, 8.5, 8.0, 5.0, 4.5, 4.0]  # 活性分子分数全部更高
    result = compute_roc_auc(labels, scores, lower_score_is_more_positive=False)
    assert result["auc"] == 1.0


def test_enrichment_factor_perfect_top_half():
    labels = [True, True, True, False, False, False]
    scores = [-9.0, -8.5, -8.0, -5.0, -4.5, -4.0]
    result = compute_enrichment_factor(labels, scores, top_fraction=0.5)
    assert result["ok"] is True
    assert result["n_actives_in_top"] == 3
    assert result["enrichment_factor"] == 2.0  # 全部活性分子集中在前50%，期望值的2倍


def test_enrichment_factor_random_order_is_close_to_one():
    """活性分子均匀散布在整个排序里，EF应该接近1(随机水平)。"""
    labels = [True, False] * 50  # 严格交替，均匀散布
    scores = list(range(100))  # 排序和标签无关
    result = compute_enrichment_factor(labels, scores, top_fraction=0.5)
    assert 0.8 < result["enrichment_factor"] < 1.2


def test_enrichment_factor_honestly_fails_without_actives():
    result = compute_enrichment_factor([False, False], [-9.0, -8.0], top_fraction=0.5)
    assert result["ok"] is False


def test_enrichment_factor_rejects_invalid_fraction():
    result = compute_enrichment_factor([True, False], [-9.0, -8.0], top_fraction=1.5)
    assert result["ok"] is False


# ------------------------------------------------------------
# compute_auc_by_similarity_to_reference：路线A的SHP2验证做完之后新增的诊断工具
# ------------------------------------------------------------
def test_similarity_grouping_detects_real_confound():
    """构造一个真实存在"混了两类机制"的场景：像苯酚的那组里分数和活性完全
    对应(真实tunnel位点逻辑)，不像苯酚的那组分数和活性完全反过来(另一种
    机制，用错了口袋)——分组后应该能看出高相似度组AUC远高于低相似度组，
    整体混在一起的AUC会被拉成很差，掩盖了真实的信号。"""
    reference = "c1ccccc1O"
    smiles = [
        "c1ccccc1O", "Cc1ccccc1O", "c1ccccc1OC", "Clc1ccccc1O",  # 像苯酚
        "CCCCCCCC", "CCN(CC)CC", "c1ccncc1", "C1CCNCC1",  # 不像苯酚
    ]
    labels = [True, True, False, False, True, True, False, False]
    scores = [-9.0, -8.5, -5.0, -5.2, -5.0, -5.1, -9.0, -8.9]
    result = compute_auc_by_similarity_to_reference(labels, scores, smiles, reference, thresholds=[0.3])
    assert result["ok"] is True
    group = result["by_threshold"][0.3]
    assert group["high_similarity_auc"]["auc"] == 1.0
    assert group["low_similarity_auc"]["auc"] == 0.0
    # 整体(不分组)的AUC应该明显比任何一组单独看都差——这正是分组前看不出真实模式的原因
    assert result["overall_auc"]["auc"] < 0.5


def test_similarity_grouping_honestly_fails_on_bad_reference():
    result = compute_auc_by_similarity_to_reference(
        [True, False], [-9.0, -8.0], ["c1ccccc1", "CCCC"], reference_smiles="not_a_valid_smiles!!!"
    )
    assert result["ok"] is False


def test_similarity_grouping_rejects_length_mismatch():
    result = compute_auc_by_similarity_to_reference(
        [True, False, True], [-9.0, -8.0], ["c1ccccc1", "CCCC"], reference_smiles="c1ccccc1O"
    )
    assert result["ok"] is False
