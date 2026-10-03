# egfr-4-drug

**没有制药/生物背景？先看 [`doc/learning-guide.md`](doc/learning-guide.md)**——
从"这个项目在解决什么问题"讲起，把EGFR/突变/激酶/共价药/降解剂/ADMET/
FEP/对接这些术语用编程类比讲清楚，读完再回来看下面的内容会顺畅很多。

奥希替尼耐药后计算驱动研发管线——三条并行路线：

| 路线 | 思路 | 目录 |
|---|---|---|
| A | SHP2/SOS1变构抑制剂，联用EGFR TKI打下游收敛点 | `routes/route_a_shp2_sos1/` |
| B | 突变选择性EGFR降解剂 | `routes/route_b_degrader/` |
| C | 四代EGFR TKI(C797S) | `routes/route_c_4th_gen_tki/` |

三条路线共用同一套计算引擎(`core/`)，区别只在于各自的配置。立项推理过程和
核实过的竞品数据见[`doc/competitive-landscape.md`](doc/competitive-landscape.md)。
**先读 [`doc/00-design-overview.md`](doc/00-design-overview.md)**——里面有
完整的三路线对比、环节→代码映射表和贯穿全仓库的工程原则。

路线C的完整科学/工程设计原文在
[`四代 EGFR (C797S) 抑制剂工业级研发管线设计.md`](<四代 EGFR (C797S) 抑制剂工业级研发管线设计.md>)，
本仓库是该设计的代码落地，后来扩展出路线A/B。

## 核心原则

凡是需要GPU集群、MD引擎、训练好的ML模型、真实专利数据这类本环境不具备的
基础设施才能产出的结果，对应函数一律诚实返回 `{"ok": False, "reason": "..."}`，
不编造数值。原因见 `core/docking.py` 的模块docstring：姊妹项目
`cancer/egfr-pipline` 的受体准备函数曾经只是 `cp receptor.pdb receptor.pdbqt`，
导致全部对接分数都是假的0.0——本仓库的每一行"真实可跑"标注都经过实际验证
（已知分子的InChIKey/分子量交叉核对、真实PDB结构编号、已上市药物分子式核对等），
每一处"诚实占位"都说明了需要什么基础设施才能真正跑起来。这条原则在三路线
扩展里延续得很彻底：路线A的SHP2结构(PDB 5EHR)、路线B的E3配体(来那度胺/
泊马度胺，分子式已核对)都经过验证，不是凭记忆编造的化学结构。

## 目录结构

```
config/
  structural_alerts.yaml       三路线共享的结构警示SMARTS规则库(已用正负对照验证，和具体靶点无关)

routes/                        路线专属配置，三条路线结构对称
  route_a_shp2_sos1/config/
    target_profile.yaml          SHP2(PDB 5EHR)为主/SOS1(PDB 6SCM)为backup，可及人群按耐药机制分型
    pipeline.yaml
    patent_landscape.yaml
    rgroup_libraries/
  route_b_degrader/config/
    target_profile.yaml          DC50/Dmax为主终点，CNS Kp,uu是一级TPP指标
    pipeline.yaml                 含三元复合物评估阶段
    patent_landscape.yaml
    rgroup_libraries/
      warhead_library.yaml        复用路线C的EGFR结合片段
      e3_ligands.yaml               来那度胺/泊马度胺衍生(已核对分子式)
      linkers.yaml                   PEG/三氮唑linker
  route_c_4th_gen_tki/config/
    target_profile.yaml          C797S基因型清单，TPP已按silevertinib真实数据校准CNS渗透target
    pipeline.yaml
    patent_landscape.yaml
    rgroup_libraries/

core/                           三路线共享的核心计算引擎
  standardize.py                   环节1 标准化 [真实可跑]
  structures.py                     环节2 结构系综 [下载/清洗真实，MD/突变占位]
  enumerate.py                       环节3 枚举 [真实可跑；固定骨架+R基团 / warhead-linker-E3两种拓扑]
  route_config.py                    路线配置路径解析(routes/<route_id>/config)
  docking.py                         环节4 L1-L2 对接 [接口正确，meeko/vina未装，占位]
  covalent.py                        环节4.3 共价对接 [攻击角判据真实，对接本身占位]
  md_stability.py / mmgbsa.py        环节4 L3 MD/MM-GBSA [占位，需GPU]
  fep.py                              环节4 L4 FEP [qc_pass逻辑真实，FEP本身占位]
  ternary_complex.py                  路线B专属：三元复合物 [hook effect逻辑真实，3D预测占位]
  selectivity.py                      环节5 选择性引擎 [热力学数学真实，依赖L4]
  admet.py                            环节6 ADMET/结构警示 [真实可跑，有精度边界声明]
  synthesis.py                        环节7 合成可及性 [SA score真实，逆合成/FTO占位]
  decide.py                           环节8 多目标决策 [真实可跑]
  feedback.py                         环节9 DMTA回传/规则治理 [真实可跑]

workflows/
  run_discovery_round.py          路线A/C共用：枚举->标准化->ADMET->决策
  run_bifunctional_round.py       路线B专属：三组分枚举->标准化->ADMET画像

tests/                           75个pytest单测，覆盖所有"真实可跑"模块+三路线配置校验+关键回归场景
db/
  schema_mysql.sql                MySQL 8.0+ schema(推荐，生产环境用这个)，三路线共用，靠route_id列区分
  schema.sql                      SQLite版，仅用于不想起MySQL server时的快速本地校验

doc/
  00-design-overview.md           总览：三路线对比 + 共享引擎架构
  01-standardization.md ... 09-feedback.md   共享引擎各模块详解(路线无关)
  competitive-landscape.md        立项决策的竞品数据核查记录(唯一事实来源)
  routes/
    route-a-shp2-sos1.md
    route-b-degrader.md
    route-c-4th-gen-tki.md         含silevertinib竞品对标的真实数据
```

## 环境搭建

```bash
cd cancer/egfr-4-drug
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

已验证可用：rdkit/pandas/numpy/pydantic/scipy/PyYAML/pytest（Python 3.9.6）。
`meeko`/`vina`（真实对接）需要Python≥3.10 + 系统Boost库，`openmm`/`openfe`/
`aizynthfinder`/`PRosettaC`（MD/FEP/逆合成/三元复合物预测）需要GPU/训练数据/
专门工具，详见 `requirements.txt` 注释。

## 快速开始

```bash
# 1. 建数据库(推荐MySQL；三条路线共用一个库，靠compounds.route_id列区分)
mysql -u root -p < db/schema_mysql.sql

# 或者不想起MySQL server时，用SQLite做快速本地校验
sqlite3 egfr4.db < db/schema.sql

# 2. 跑全部测试
.venv/bin/python -m pytest tests/ -v

# 3. 路线C端到端工作流(默认route_id)：枚举324个候选分子 -> 标准化 ->
#    ADMET/结构警示硬门 -> Pareto前沿 -> 批次选择，几秒内跑完
.venv/bin/python -m workflows.run_discovery_round

# 4. 路线A端到端工作流(同一个workflow，换route_id)
.venv/bin/python -c "from workflows.run_discovery_round import run; run(route_id='route_a_shp2_sos1', batch_size=10)"

# 5. 路线B端到端工作流(三组分枚举+ADMET画像；决策阶段诚实缺位，见该文件docstring)
.venv/bin/python -m workflows.run_bifunctional_round

# 6. 单独跑某个模块的演示(每个core/*.py都有可独立运行的__main__演示)
.venv/bin/python -m core.standardize
.venv/bin/python -m core.admet
.venv/bin/python -m core.decide
.venv/bin/python -m core.selectivity
.venv/bin/python -m core.feedback
.venv/bin/python -m core.ternary_complex
.venv/bin/python -m core.route_config
```

## 文档索引

| 主题 | 文档 |
|---|---|
| **新手入门(无制药背景)** | [doc/learning-guide.md](doc/learning-guide.md) |
| 总览(三路线对比) | [doc/00-design-overview.md](doc/00-design-overview.md) |
| 竞品数据核查记录 | [doc/competitive-landscape.md](doc/competitive-landscape.md) |
| 路线A：SHP2/SOS1 | [doc/routes/route-a-shp2-sos1.md](doc/routes/route-a-shp2-sos1.md) |
| 路线B：EGFR降解剂 | [doc/routes/route-b-degrader.md](doc/routes/route-b-degrader.md) |
| 路线C：四代EGFR TKI | [doc/routes/route-c-4th-gen-tki.md](doc/routes/route-c-4th-gen-tki.md) |
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
   ——已修复并加回归测试（见 `doc/03-enumeration.md`）。路线B新增的
   warhead-linker-E3三组分拼接复用了同一套连通性校验逻辑。
4. **选择性引擎符号约定写反**（`core/selectivity.py`文档字符串）——已修复并
   用无歧义场景钉死方向（见 `doc/05-selectivity.md`）。
5. **`standardize_batch()` 在"全成功"批次下会崩溃**——DataFrame缺列导致
   KeyError，端到端集成测试才暴露，单模块demo测不出来（见
   `doc/01-standardization.md`）。
6. **`target_profile.yaml` 的 `decision_engine.pareto_objectives` 引用了
   tpp表里不存在的指标**（`cns_mpo`/`synthesis_n_steps`），且`core/decide.py`
   从未真正读取过这个配置——已补全tpp表并实现对接函数（见
   `doc/08-decision.md`）。
7. **路线C的`chemistry_route`注释把precedent搞错了**：原文说"主攻可逆正构
   的理由是机制已被silevertinib验证"，但核查后确认silevertinib实际是共价药，
   真正的可逆正构precedent是tigozertinib(BLU-945)——已修正，同时把
   `chemistry_route`枚举从字母(A/B/C/D)改成描述性slug，避免和本次新增的
   "路线A/B/C"三路线框架产生命名混淆（见
   `doc/routes/route-c-4th-gen-tki.md`）。

## 数据库选型：MySQL/Redis/ES是否需要接入

当前管线的数据量级（一轮几十到几百个化合物，三条路线累计到项目生命周期内
也就几千到几万行）用MySQL单机绰绰有余，`db/schema_mysql.sql` 已经是可以
直接生产使用的版本。Redis和Elasticsearch**现阶段不需要接入**：

- **Redis**：管线里没有高频重复计算或分布式任务队列的场景——`core/admet.py`/
  `core/synthesis.py` 这些纯函数按InChIKey缓存的收益在当前规模可以忽略。等以后
  真的接入Vina/MD/FEP、需要把L1-L4打分漏斗分发给多个worker并行跑时(三条路线
  一起跑尤其明显)，Redis做Celery/RQ的broker才有意义。
- **Elasticsearch**：我们的核心查询是结构化聚合(Pareto前沿、按scaffold分组的
  MUE)，SQL足够。ES真正有价值的场景是化合物库涨到10万+后用dense_vector/kNN
  做分子指纹相似度搜索——不是现在这个阶段的需求。
