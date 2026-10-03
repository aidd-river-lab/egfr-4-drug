# 环节9：DMTA闭环、主动学习与规则治理

代码：[`core/feedback.py`](../core/feedback.py) · **状态：真实可跑**

## 为什么这一环节存在

前面所有环节产出的都是**预测**。环节9是整条管线唯一能验证这些预测对不对的
地方——没有这一环，管线只是一套自洽但可能系统性跑偏的假设。原设计文档
"环节9.2"规定了每轮DMTA(Design-Make-Test-Analyze)循环必做的四件事，本模块
实现了其中可以用纯代码表达的部分。

## 湿实验回传Schema（环节9.1）

`WetLabResult` 用pydantic实现（不是自由格式Excel），赋值时就强校验类型：

```python
from core.feedback import WetLabResult
from datetime import date

WetLabResult(
    compound_id="INT-0412", batch_id="B001",
    smiles_as_made="...",           # 实际做出来的结构，不能省
    purity_pct=97.5, chirality_confirmed=True,
    genotype="del19_C797S", assay_format="cellular",
    readout="pIC50", value=8.1, unit="pIC50", n_replicates=3,
    assay_date=date(2026, 10, 1), operator="...", protocol_version="v1",
)
```

字段分六组，逐字照抄原设计文档"环节9.1"的表格：

- **身份**：compound_id/batch_id/**smiles_as_made**/purity_pct/chirality_confirmed
- **活性**：genotype/assay_format/atp_conc_um/readout/value/unit/censoring/
  n_replicates/log_sd/assay_date/plate_id
- **选择性**：wt_value/fold_selectivity/kinase_panel_file
- **ADME**：sol_ph74/hlm_clint/rlm_clint/ppb_pct/caco2_papp/mdr1_er/herg_ic50_um/
  gsh_adduct_detected
- **体内**：species/dose_mg_per_kg/route/auc/cmax/t_half_hours/f_pct/kp_uu
- **元数据**：operator/protocol_version/qc_status/comments

**`smiles_as_made` 为什么不能和"设计出来的结构"合并成一个字段**：合成出来的
实际结构和设计稿不一致是家常便饭(立体化学翻转、副产物、不完全反应产物被
误认为主产物)，如果拿实测活性去对设计结构做训练/复盘，整个预测模型会被
系统性污染，而且这种污染很难事后发现——表现上看就是"模型莫名其妙预测不准"，
但根因是数据输入错了，不是模型本身的问题。

## 预测vs实测复盘（环节9.2第1条）

```python
from core.feedback import compare_predicted_vs_actual

compare_predicted_vs_actual(df, predicted_col="predicted_pic50", actual_col="actual_pic50", scaffold_col="scaffold")
# {"spearman_rho": 0.89, "mue": 0.35, "per_scaffold": {"series_A": {"mue": 0.32}, "series_B": {"mue": 0.38}}}
```

整体Spearman ρ和MUE（真实用`scipy.stats.spearmanr`计算），以及按化学型分组的
误差——原设计文档强调"整体ρ好但某个子系列系统性偏差"是最有价值的信号之一，
因为它经常指向"这个子系列有某种结构特征(比如特定的构象限制)是模型没学到的"，
而不是模型整体失效。

## 假阳性归因（环节9.2第2条）：结构化标签，不是自由文本

```python
FalsePositiveAttribution = Literal[
    "pose_error", "protonation_state_error", "conformer_selection_error",
    "desolvation_underestimated", "cell_permeability_issue",
    "metabolic_instability", "compound_degradation_or_purity",
]
```

原设计文档强调"标签要进数据库，季度统计"——如果用自由文本记录归因原因，
没法做季度统计也没法训练"哪类归因在哪类骨架上更常见"这种二阶分析。落库表：
`db/schema.sql` 的 `false_positive_attributions`，`attribution_label`列用
`CHECK`约束强制只能是这7个值之一。

## 规则治理（环节9.3）：precision追踪决定规则能不能升级成硬拒绝

```python
from core.feedback import rule_precision_tracking

rule_precision_tracking(rule_hits_df, true_label_col="actually_bad", rule_hit_col="rule_hit")
# {"precision": 0.857, "n_independent_compounds": 7, "eligible_for_reject_mode": True, ...}
```

一条结构警示规则要从"警告(OPTIMIZE_NEEDED)"模式升级成"硬拒绝(REJECT)"模式，
必须同时满足：precision > 0.8 **且**在至少3个独立化合物上验证过——这是
原设计文档"环节9.3"的硬要求，防止一条规则因为偶然在1-2个化合物上命中就被
过早固化成拒绝规则。落库表：`db/schema.sql` 的 `rule_precision_history`，
`current_mode`列记录规则当前实际处于哪个模式。

## 和环节8的闭环

`core/decide.py::select_batch()` 的 `rejected_pool` 参数（被规则/低分拒绝的
候选池）配合本模块的 `rule_precision_tracking()`，形成完整的"假阴性复活"
闭环：规则拒绝的分子一部分会被 `select_batch()` 的 `control` 档随机抽样
重新送去合成和测试，测试结果反过来喂给 `rule_precision_tracking()` 统计这条
规则的真实precision，决定要不要继续信任它、或者要不要把它从警告升级成拒绝。

## 前瞻验证（环节9.2第4条）

`db/schema.sql` 的 `prospective_validation_rounds` 表记录每一轮模型重训后的
整体预测质量（Spearman ρ/MUE随时间的序列），用来判断模型是在变好还是在退化。
这部分目前只有落库schema，没有单独的分析函数——它本质是对
`compare_predicted_vs_actual()` 的结果按轮次做时间序列记录，不需要额外的
计算逻辑，接入时直接把每轮调用结果写进这张表即可。
