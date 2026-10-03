# 环节8：多目标决策引擎

代码：[`core/decide.py`](../core/decide.py) · 测试：[`tests/test_decide.py`](../tests/test_decide.py) · **状态：真实可跑**

## 为什么不用加权求和

原设计文档"环节8.1"明确反对加权求和："高分可以靠一项极端值补偿其他项的缺陷"
——一个hERG风险很高但细胞活性极好的分子，加权求和打分可能排名很靠前，但这样
的分子根本不应该被合成。本模块实现两条推荐路线，都不允许这种"拆东墙补西墙"。

## 1. Pareto前沿（非支配排序）

```python
from core.decide import pareto_front

pareto_front(df, {"cell_ic50_nm": "lower_is_better", "wt_selectivity_fold": "higher_is_better", "cns_mpo": "higher_is_better"})
# 返回布尔Series：True表示该行是非支配解(不存在另一行在所有目标上都不比它差、且至少一项严格更好)
```

纯numpy广播实现，不依赖任何第三方Pareto库。

## 2. Desirability score（几何平均）

```python
from core.decide import DesirabilitySpec, desirability_score

specs = [
    DesirabilitySpec("cell_ic50_nm", lambda x: 1.0 if x <= 20 else max(0.0, 1 - (x - 20) / 80)),
    DesirabilitySpec("wt_selectivity_fold", lambda x: min(1.0, x / 100)),
    DesirabilitySpec("cns_mpo", lambda x: min(1.0, x / 4.5)),
]
desirability_score(df, specs)
```

每个指标先算0-1的desirability分数，取加权**几何平均**作为总分。几何平均的
关键性质：任何一项desirability=0会让总分直接变成0——这是有意的设计，不是bug。
`tests/test_decide.py::test_desirability_is_not_a_weighted_sum` 专门验证了
这一点：一项1.0分和一项0.0分的几何平均必须是0，不能像加权和那样错误地算出0.5。

## 3. 批次选择：利用/探索/假设检验/对照（环节8.3）

```python
from core.decide import select_batch

batch = select_batch(
    candidates, desirability_col="desirability", batch_size=40,
    uncertainty_col="uncertainty",              # 没有就不传，explore档会诚实留空
    hypothesis_test_indices=[5, 12, 31],        # 药化指定，不传就留空
    rejected_pool=rejected_df,                  # 假阴性复活抽样来源
)
batch.exploit / batch.explore / batch.hypothesis_test / batch.control  # 四个DataFrame
batch.combined  # 合并后带batch_role列
```

默认配额50/25/15/10（`config/pipeline.yaml` 的 `batch_composition`）：

- **exploit**：按desirability取top-N，用期望效用最大化。
- **explore**：按不确定度取top-N，用UCB类采集函数——**没有不确定度列就跳过，
  不编一个假数据去填**（真实的不确定度要来自ML ensemble方差或FEP误差估计）。
- **hypothesis_test**：这一档本质是药化的化学判断（原设计文档例子："测试797
  入口区砜基是否必要"这种成对假设），**不是算法能自动生成的**——只接受调用方
  显式传入的index列表。不传就留空，**配额不会被挪用给其他档位去凑数**，
  免得掩盖"这一轮没做假设检验"这个事实。`tests/test_decide.py::test_select_batch_does_not_reallocate_empty_hypothesis_quota`
  验证了总选出数量会诚实地少于batch_size，而不是悄悄把配额分给别的档位。
- **control**：优先从 `rejected_pool`（被规则/低分拒绝的候选池）随机抽样，
  用来检验规则有没有误杀好分子——这是原设计文档"环节8.3/9.3"反复强调的
  "假阴性复活机制"的具体实现，和 `core/feedback.py` 的
  `rule_precision_tracking()` 配合，形成"规则拒绝→部分样本仍进入control组
  测试→事后统计precision→决定规则是否可信"的闭环。

## 和 `target_profile.yaml` 的对接

`target_profile.yaml` 里有一条注释："决策引擎读取这里的字段名来对齐TPP...
不要改字段名除非同步改TPP_FIELD_MAP"——这条注释曾经只是声明，没有对应代码，
配置和代码实际上没接上。现在通过 `load_decision_config_from_target_profile()`
把这个对接点补上了：

```python
from core.decide import load_decision_config_from_target_profile

load_decision_config_from_target_profile()
# {
#   "pareto_objectives": {"cell_ic50_del19_c797s_nm": "lower_is_better",
#                         "wt_selectivity_fold": "higher_is_better",
#                         "cns_mpo": "higher_is_better",
#                         "synthesis_n_steps": "lower_is_better"},
#   "missing_from_tpp_table": [],
#   "skipped_non_monotonic_metrics": [],
#   "batch_quota": {"exploit_pct": 50, "explore_pct": 25, "hypothesis_test_pct": 15, "control_pct": 10},
# }
```

实现过程中发现并修复了一个真实的配置bug：`decision_engine.pareto_objectives`
列出的 `cns_mpo` 和 `synthesis_n_steps` 两个指标，最初根本没有在 `tpp` 表里
登记对应的行——这正是 `target_profile.yaml` 文件顶部自己警告的"不允许有打分
维度在这里找不到对应项，那是'自欺'的信号"。修复方式是给这两个指标各补了一行
完整的TPP定义(minimum/target/direction/unit/computed_by)。函数会把"指标在
tpp表里完全缺失"(`missing_from_tpp_table`，真正的配置bug，应该被修掉)和
"指标存在但方向是`in_range`不能直接喂给Pareto"(`skipped_non_monotonic_metrics`，
比如`clogd74`，这是预期行为)分开报告，不要混为一谈。

## 端到端用法

`workflows/run_discovery_round.py` 串联了本环节和环节1/3/6/7，是目前唯一一个
完整跑通的端到端工作流，见 [README](../README.md) 的快速开始部分。
