# 设计总览

本文档是整个仓库文档集的入口。完整的科学/工程设计原文在
[`四代 EGFR (C797S) 抑制剂工业级研发管线设计.md`](<四代 EGFR (C797S) 抑制剂工业级研发管线设计.md>)
（以下简称"原设计文档"），本文件只负责说清楚：这份设计被落地成了哪些代码、
代码在哪、以及每一部分现在是"真实可跑"还是"诚实占位"。

## 项目背景

奥希替尼（osimertinib，三代EGFR抑制剂）获得性耐药后最常见的机制之一是EGFR激酶域
新增 C797S 突变——它恰好发生在奥希替尼共价成键的那个半胱氨酸上，让第三代药物
彻底失效。本项目的目标是设计能克服 C797S 耐药、同时保留对野生型EGFR选择性
（避免皮疹/腹泻等剂量限制性毒性）的第四代抑制剂。

详见 `config/target_profile.yaml` 的"0.1 基因型清单"：主攻 del19/C797S 和
L858R/C797S（无T790M，这是目前最大也最难的病人群），次要覆盖两个三突变体，
反靶点包括 EGFR 野生型和 7 个脱靶激酶。

## 管线结构：10个环节

原设计文档把管线分成环节0-9，每个环节在本仓库里对应的落地状态：

| 环节 | 名称 | 代码位置 | 状态 | 详细文档 |
|---|---|---|---|---|
| 0 | 目标画像(TPP) | `config/target_profile.yaml` | 配置，真实 | 本文件 |
| 1 | 分子标准化 | `core/standardize.py` | **真实可跑** | [01](01-standardization.md) |
| 2 | 结构系综准备 | `core/structures.py` | 部分真实+占位 | [02](02-structure-ensemble.md) |
| 3 | 类似物枚举 | `core/enumerate.py` | **真实可跑** | [03](03-enumeration.md) |
| 4 | 多级打分漏斗 L0-L4 | `core/docking.py` `core/mmgbsa.py` `core/md_stability.py` `core/fep.py` `core/covalent.py` | 部分真实+占位 | [04](04-funnel-docking.md) |
| 5 | 选择性引擎 | `core/selectivity.py` | 真实数学+占位计算 | [05](05-selectivity.md) |
| 6 | ADMET/结构警示 | `core/admet.py` | **真实可跑**(有精度边界声明) | [06](06-admet.md) |
| 7 | 合成可及性 | `core/synthesis.py` | 部分真实(SA score)+占位 | [07](07-synthesis.md) |
| 8 | 多目标决策 | `core/decide.py` | **真实可跑** | [08](08-decision.md) |
| 9 | DMTA回传与规则治理 | `core/feedback.py` | **真实可跑**(schema+分析逻辑) | [09](09-feedback.md) |

## 核心工程原则（贯穿所有模块，不是口号）

1. **不编造数字**。凡是需要GPU集群、MD引擎、训练好的ML模型、真实专利数据库
   这类本环境不具备的基础设施才能产出的结果，对应函数一律返回
   `{"ok": False, "reason": "..."}`，把缺什么说清楚，不返回一个"看起来合理"
   的假数值。这个原则的由来：姊妹项目 `cancer/egfr-pipline` 的受体PDBQT
   准备函数实际上是 `cp receptor.pdb receptor.pdbqt`（只是换后缀名，完全
   没有做AutoDock需要的原子类型/电荷计算），导致那个项目产出的全部对接
   分数都是 0.0 kcal/mol——一个看起来能跑、实际全错的典型反面教材。
   `core/docking.py` 的模块docstring里专门记录了这个教训。

2. **成药性是乘法不是加法**。`core/decide.py` 的多目标打分用几何平均
   (desirability score)而不是加权和，任何一项指标为0会让总分归零，
   不允许靠某一项极端值掩盖其他项的致命缺陷。

3. **规则要有豁免清单，且豁免清单要有验证过的正负对照**。
   `config/structural_alerts.yaml` 里每条SMARTS规则都用真实分子验证过
   （奥希替尼、吉非替尼作为正常分子不应该被误报），验证脚本见
   `tests/test_structural_alerts.py`。

4. **硬门槛要慎用，弱信号不升级成拒绝规则，除非precision经过验证**。
   `core/feedback.py` 的 `rule_precision_tracking()` 实现了这条治理机制：
   一条规则要从"警告"升级成"拒绝"，必须precision>0.8且在≥3个独立化合物上
   验证过。

5. **字段schema要和设计文档的表格逐字对齐**，下游代码按字段名读取，
   不允许每个模块自己发明一套命名。`db/schema.sql` 的每张表、
   `core/*.py` 里的每个dataclass/pydantic模型都逐字照抄了原设计文档对应
   小节的字段表。

## 环境与依赖

见 [`requirements.txt`](../requirements.txt)。核心依赖（RDKit/pandas/numpy/
pydantic/scipy/PyYAML/pytest）在本仓库 `.venv`（Python 3.9.6）下已验证可用。
`meeko`/`vina`（对接需要）和 `openmm`/`openfe`/`aizynthfinder`（MD/FEP/逆合成
需要）未安装，原因和升级路径见该文件注释。

## 快速开始

```bash
# 建库
sqlite3 egfr4.db < db/schema.sql

# 跑测试(47+个单测，覆盖所有真实可跑的模块)
.venv/bin/python -m pytest tests/ -v

# 跑端到端工作流(环节1+3+6+7+8的完整一轮，纯RDKit，不需要GPU)
.venv/bin/python -m workflows.run_discovery_round

# 单独跑某个模块的演示
.venv/bin/python -m core.standardize
.venv/bin/python -m core.admet
.venv/bin/python -m core.decide
```

## 已知的、故意保留的局限

- `core/admet.py` 的 `clogd_approx` 直接复用 `clogp`(中性分子分配系数)，
  不是真正按pH7.4电离态算的cLogD——需要pKa预测器才能做对，模块顶部有
  完整声明。
- `config/rgroup_libraries/scaffolds.yaml` 里只有一个示例性骨架
  (`demo_aminopyrimidine_biphenyl`)，不是真实候选骨架，仅用来演示数据结构。
- `config/patent_landscape.yaml` 是占位数据(`status: placeholder_not_real`)，
  不能用于真实FTO决策。
- 环节2/4(L1往后)/5/9 的真实计算都需要本环境不具备的基础设施(PDB结构库+
  MD软件+GPU集群、湿实验数据)，当前是正确接口+诚实占位，接入方式见各自
  的分文档。
