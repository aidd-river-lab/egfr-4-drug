# egfr-4-drug

EGFR C797S 第四代抑制剂（奥希替尼耐药后）计算驱动研发管线。

完整科学/工程设计见
[`四代 EGFR (C797S) 抑制剂工业级研发管线设计.md`](<四代 EGFR (C797S) 抑制剂工业级研发管线设计.md>)，
本仓库是该设计的代码落地。**先读 [`doc/00-design-overview.md`](doc/00-design-overview.md)**
——里面有完整的环节→代码映射表和贯穿全仓库的工程原则。

## 核心原则

凡是需要GPU集群、MD引擎、训练好的ML模型、真实专利数据这类本环境不具备的
基础设施才能产出的结果，对应函数一律诚实返回 `{"ok": False, "reason": "..."}`，
不编造数值。原因见 `core/docking.py` 的模块docstring：姊妹项目
`cancer/egfr-pipline` 的受体准备函数曾经只是 `cp receptor.pdb receptor.pdbqt`，
导致全部对接分数都是假的0.0——本仓库的每一行"真实可跑"标注都经过实际验证
（已知分子的InChIKey/分子量交叉核对、合成数据的边界场景测试等），每一处
"诚实占位"都说明了需要什么基础设施才能真正跑起来。

## 目录结构

```
config/                     环节0：唯一真相来源(TPP/基因型/化学路线/规则/批次配额)
  target_profile.yaml          目标画像(TPP表、基因型清单、决策引擎配置)
  structural_alerts.yaml       结构警示SMARTS规则库(已用正负对照验证)
  patent_landscape.yaml        专利占位数据(明确标注非真实，勿用于FTO决策)
  pipeline.yaml                L0-L4打分漏斗配置、共价路线、批次配额
  rgroup_libraries/
    scaffolds.yaml              骨架模板(示例性)
    pocket_regions.yaml         按口袋区域组织的R基团库

core/                        各环节的核心计算逻辑
  standardize.py               环节1 标准化 [真实可跑]
  structures.py                环节2 结构系综 [下载/清洗真实，MD/突变占位]
  enumerate.py                  环节3 枚举 [真实可跑]
  docking.py                    环节4 L1-L2 对接 [接口正确，meeko/vina未装，占位]
  covalent.py                   环节4.3 共价对接 [攻击角判据真实，对接本身占位]
  md_stability.py / mmgbsa.py   环节4 L3 MD/MM-GBSA [占位，需GPU]
  fep.py                        环节4 L4 FEP [qc_pass逻辑真实，FEP本身占位]
  selectivity.py                环节5 选择性引擎 [热力学数学真实，依赖L4]
  admet.py                      环节6 ADMET/结构警示 [真实可跑，有精度边界声明]
  synthesis.py                  环节7 合成可及性 [SA score真实，逆合成/FTO占位]
  decide.py                     环节8 多目标决策 [真实可跑]
  feedback.py                   环节9 DMTA回传/规则治理 [真实可跑]

workflows/
  run_discovery_round.py       端到端工作流：串联环节1+3+6+7+8(纯RDKit，不需要GPU)

tests/                        47+个pytest单测，覆盖所有"真实可跑"模块+关键回归场景
db/
  schema.sql                   SQLite schema，17张表，字段名对齐各core/*.py的dataclass

doc/                          本文档集，每个环节一篇
  00-design-overview.md ... 09-feedback.md
```

## 环境搭建

```bash
cd cancer/egfr-4-drug
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

已验证可用：rdkit/pandas/numpy/pydantic/scipy/PyYAML/pytest（Python 3.9.6）。
`meeko`/`vina`（真实对接）需要Python≥3.10 + 系统Boost库，`openmm`/`openfe`/
`aizynthfinder`（MD/FEP/逆合成）需要GPU/训练数据，详见 `requirements.txt` 注释。

## 快速开始

```bash
# 1. 建数据库
sqlite3 egfr4.db < db/schema.sql

# 2. 跑全部测试
.venv/bin/python -m pytest tests/ -v

# 3. 跑端到端工作流：枚举324个候选分子 -> 标准化 -> ADMET/结构警示硬门 ->
#    Pareto前沿 -> 批次选择，全程纯RDKit，几秒内跑完
.venv/bin/python -m workflows.run_discovery_round

# 4. 单独跑某个模块的演示(每个core/*.py都有可独立运行的__main__演示)
.venv/bin/python -m core.standardize
.venv/bin/python -m core.admet
.venv/bin/python -m core.decide
.venv/bin/python -m core.selectivity
.venv/bin/python -m core.feedback
```

## 文档索引

| 环节 | 文档 |
|---|---|
| 总览 | [doc/00-design-overview.md](doc/00-design-overview.md) |
| 1 标准化 | [doc/01-standardization.md](doc/01-standardization.md) |
| 2 结构系综 | [doc/02-structure-ensemble.md](doc/02-structure-ensemble.md) |
| 3 枚举 | [doc/03-enumeration.md](doc/03-enumeration.md) |
| 4 打分漏斗 | [doc/04-funnel-docking.md](doc/04-funnel-docking.md) |
| 5 选择性 | [doc/05-selectivity.md](doc/05-selectivity.md) |
| 6 ADMET | [doc/06-admet.md](doc/06-admet.md) |
| 7 合成可及性 | [doc/07-synthesis.md](doc/07-synthesis.md) |
| 8 多目标决策 | [doc/08-decision.md](doc/08-decision.md) |
| 9 DMTA回传 | [doc/09-feedback.md](doc/09-feedback.md) |

## 开发中发现并修复的真实bug（如果要继续开发，先读这些）

1. **奥希替尼SMILES错误**导致姊妹项目全部相似度算成0.0——根因是参考分子
   SMILES无法通过RDKit kekulize，已找到并验证正确的SMILES（见
   `doc/01-standardization.md`）。
2. **结构警示规则误报奥希替尼**（`aniline_unmasked`/`diaminobenzene_1_4/1_2`）
   ——已修复并用吉非替尼交叉验证（见 `doc/06-admet.md`）。
3. **molzip同位素不匹配导致断裂分子静默通过SanitizeMol**（`core/enumerate.py`）
   ——已修复并加回归测试（见 `doc/03-enumeration.md`）。
4. **选择性引擎符号约定写反**（`core/selectivity.py`文档字符串）——已修复并
   用无歧义场景钉死方向（见 `doc/05-selectivity.md`）。
5. **`standardize_batch()` 在"全成功"批次下会崩溃**——DataFrame缺列导致
   KeyError，端到端集成测试才暴露，单模块demo测不出来（见
   `doc/01-standardization.md`）。
6. **`target_profile.yaml` 的 `decision_engine.pareto_objectives` 引用了
   tpp表里不存在的指标**（`cns_mpo`/`synthesis_n_steps`），且`core/decide.py`
   从未真正读取过这个配置——已补全tpp表并实现对接函数（见
   `doc/08-decision.md`）。
