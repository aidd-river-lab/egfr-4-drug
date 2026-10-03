# 环节1：分子标准化

代码：[`core/standardize.py`](../core/standardize.py) · 测试：[`tests/test_standardize.py`](../tests/test_standardize.py) · **状态：真实可跑**

## 为什么需要这一步

同一个分子可以写成很多种不同但化学上等价的SMILES：带不带盐、互变异构体写法、
电荷表示方式都可能不同。如果不统一，后续所有基于SMILES字符串或简单指纹的去重、
聚合、建模都会把同一个化合物误判成不同的化合物（或反过来，把不同化合物误判成
重复）。环节1把这件事一次做对，下游所有环节都依赖它产出的 `inchikey`/`inchikey14`
做身份识别。

## 流程（严格按顺序，顺序本身是设计的一部分）

1. `Chem.MolFromSmiles(sanitize=True)` 解析——失败的分子进隔离表，不静默丢弃。
2. `LargestFragmentChooser` 去盐，只保留最大的有机片段。
3. `Uncharger` 中和形式电荷——注意这一步**不**固定质子化态，只是把能中和的电荷
   中和掉，真实的pH依赖质子化态需要专门的pKa预测器（见 `core/admet.py` 的
   `pka_proxy` 声明）。
4. `TautomerEnumerator.Canonicalize` 规范化互变异构体——这一步比去盐更容易被
   忽略，但吡唑/酰胺这类互变异构体如果不规范化，同一个化合物会产出两个不同的
   InChIKey，直接污染后续的去重和SAR聚合。
5. 立体化学全程不触碰——`LargestFragmentChooser`/`Uncharger`/
   `TautomerEnumerator` 三步都不会清除手性标记，标准化后用
   `Chem.AssignStereochemistry` 做一次合法性复查。
6. 生成两套主键：完整 `InChIKey`（精确去重，含立体信息）和 `InChIKey14`
   （只取前14位，去立体，用于骨架级聚合——比如"这些对映异构体是不是来自同一个
   设计思路"这种问题）。

## 接口

```python
from core.standardize import standardize_batch

result = standardize_batch(["CCO", "not-a-smiles", "c1ccccc1.[Na+]"])
result.curated      # DataFrame：标准化成功的行
result.quarantine   # DataFrame：失败的行 + 失败原因(fail_stage/fail_reason)
result.success_rate # 成功率
```

`standardize_one()` 是单分子版本，返回一个 dict，字段和上面 DataFrame 的列一致。

## 已验证的正确性

- 奥希替尼的 InChIKey 计算结果为 `DUYJMQONPNNFPI-UHFFFAOYSA-N`，与公开数据库一致
  （`tests/test_standardize.py::test_osimertinib_inchikey_matches_known_value`）。
- `c1ccccc1.[Na+].[Cl-]`（苯+氯化钠盐）正确地只剩下苯。
- 非法SMILES正确进入隔离表而不是被吞掉。

## 一个端到端集成才暴露出来的bug（已修复）

`standardize_batch()` 早期实现在"全部成功"或"全部失败"这种批次里会直接崩溃：
`pd.DataFrame(rows)` 由一组 dict 构造，如果没有任何一行带 `fail_stage`/
`fail_reason` 这两个key（全部成功时就是这样），pandas 根本不会生成这两列，
后面按列名取子集 `df[~df["ok"]][["raw_smiles","fail_stage","fail_reason"]]`
就会抛 `KeyError`。

单独跑模块自带的demo（故意混合了成功和失败样本）测不出这个问题，是在接入
`workflows/run_discovery_round.py` 做真实的324分子批量标准化（这一批全部标准化
成功，quarantine为空）时才暴露出来的。修复方式是在构造DataFrame后显式补全这两列
（`core/standardize.py` 里有详细注释），并补了三个回归测试锁定"全成功"、
"全失败"、"混合"三种情况都不崩溃。这是一个很好的例子说明：**模块自己的demo
通过不代表集成没问题**，必须跑真实的端到端数据流才能发现这类边界条件bug。
