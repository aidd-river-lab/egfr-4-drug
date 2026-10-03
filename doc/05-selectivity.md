# 环节5：选择性引擎

代码：[`core/selectivity.py`](../core/selectivity.py) · 测试：[`tests/test_selectivity.py`](../tests/test_selectivity.py) · **状态：热力学数学真实可跑，FEP计算本身依赖环节4 L4(占位)**

原设计文档原话："整条管线最该投入算力的地方"——因为WT选择性(`wt_selectivity_fold`)
直接决定皮疹/腹泻这类剂量限制性毒性的严重程度，是TPP表里权重很高的一项
（`target_profile.yaml` 的 `wt_selectivity_fold` 目标是100倍）。

## 方法：炼金术双微扰（Alchemical Double Perturbation）

不是分别对突变体和WT各算一次绝对结合自由能再相减（误差会叠加到无法使用的
量级），而是在突变体和WT**两套结构**上跑**同一组配体微扰**(L1→L2的化学改造)，
然后用两次相对结合自由能(RBFE)相减得到ΔΔΔG：

```
RBFE(L1→L2) = ΔG_L2 - ΔG_L1        # FEP的标准输出，"更负=结合更有利"
ΔΔΔG = RBFE_突变体(L1→L2) - RBFE_WT(L1→L2)
```

```python
from core.selectivity import compute_dddg_from_fep_maps
from core.fep import FEPMapResult, FEPPerturbation

mutant_map = FEPMapResult(genotype="del19_C797S", perturbations=[
    FEPPerturbation("compound_A", "compound_B", predicted_ddg_kcal_mol=-1.5, uncertainty_kcal_mol=0.3),
])
wt_map = FEPMapResult(genotype="EGFR_WT", perturbations=[
    FEPPerturbation("compound_A", "compound_B", predicted_ddg_kcal_mol=0.2, uncertainty_kcal_mol=0.25),
])
compute_dddg_from_fep_maps("compound_A", "compound_B", "del19_C797S", mutant_map, wt_map)
# {"dddg_kcal_mol": -1.7, "fold_selectivity_pred": 17.6, ...}
```

要求两张FEP map里必须能找到完全一样的微扰对，否则无法相减——这是方法学上的
硬要求（两次独立的绝对自由能没有意义，必须是同一个化学变换），不是实现偷懒；
找不到时诚实返回 `ok=False`。

## 符号约定（容易搞反，本模块踩过一次这个坑）

`ΔΔΔG < 0` 表示"L1→L2这个改造带来的结合增益，在突变体里比在WT里更大"，
也就是L2相对L1对突变体更有选择性偏好，`fold_selectivity_pred > 1`。
`ΔΔΔG > 0` 则相反：改造让分子更偏向WT，选择性变差，`fold < 1`。

**这个方向在第一版文档字符串里写反过一次**：误写成"ΔΔΔG>0表示对突变体更有利"，
是在对照实际验证过的demo输出（`dddg=-1.7 → fold=17.6`，明显是"偏向突变体"）
时才发现矛盾并改正的。为了不让这类符号错误再次悄悄溜回来，
`tests/test_selectivity.py` 专门用一个无歧义的合成场景钉死方向：

```python
def test_unambiguous_synthetic_scenario_mutant_strongly_favored():
    # 同一个改造：突变体上结合大幅变好(RBFE=-3.0)，WT上完全不变(RBFE=0.0)
    # 这种改造显然应该让fold >> 1，如果代码算出fold < 1，说明符号又被绕反了
    ...
    assert result["dddg_kcal_mol"] < 0
    assert result["fold_selectivity_pred"] > 100
```

## 热力学转换

基于 ΔG = -RT ln(K) 的玻尔兹曼关系（T=298.15K）：

```python
def dddg_to_fold_selectivity(dddg_kcal_mol, temperature_k=298.15):
    return math.exp(-dddg_kcal_mol / (R_KCAL_PER_MOL_K * temperature_k))
```

两次独立RBFE相减时的不确定度用标准误差传播公式（平方和开根号）：

```python
def propagate_subtraction_uncertainty(u1, u2):
    return math.sqrt(u1**2 + u2**2)
```

这两个函数都是纯数学，不需要真的跑过FEP就能完整验证正确性——这是本模块
"可以独立于环节4 L4先测试、先用"的部分。

## 数据模型

`SelectivityResult` 字段完全照抄原设计文档"环节5.4"给出的接口形状：

```python
@dataclass
class SelectivityResult:
    ligand_id: str
    ddg_target: dict[str, float]         # 每个目标基因型的ΔG
    ddg_wt: float | None
    dddg: dict[str, float]               # 选择性差值，排序用主指标
    fold_selectivity_pred: dict[str, float]
    uncertainty: dict[str, float]
    method: Literal["fep", "mmgbsa", "ml", "docking"]
    structures_used: list[str]
```

落库表：`db/schema.sql` 的 `selectivity_results`（每个化合物×基因型×方法一行）。

## 依赖关系

本模块的真实计算输出依赖环节4 L4(FEP)产出的 `FEPMapResult`，而L4本身在本环境
是诚实占位（见 [04-funnel-docking.md](04-funnel-docking.md)）。也就是说，
`core/selectivity.py` 的数学逻辑已经完整实现并测试通过，一旦环节4 L4在有GPU的
环境跑出真实FEP map，选择性计算可以直接接入，不需要改这个模块的任何代码。
