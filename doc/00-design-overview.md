# 设计总览

本文档是整个仓库文档集的入口。**推荐的阅读起点，按顺序：**

1. **[learning-guide.md](learning-guide.md)**——零生化背景的术语速成，
   用编程类比把EGFR、突变、共价药、降解剂、ADMET、对接/FEP这些术语讲清楚。
2. **[16-visual-dashboard.md](16-visual-dashboard.md)**——**一张图看懂
   全局**：问题背景、三路线战略、管线流程、三路线真实进度、已验证的结果、
   和工业级的差距、接下来怎么走，全部配图，不需要先读完其他文档。
3. 以上两篇看完之后，按兴趣挑：
   - 想知道代码具体怎么跑、黑盒打分函数背后是什么数学——
     [10-pipeline-walkthrough.md](10-pipeline-walkthrough.md)
   - 想跟着两个真实分子(RTC-0000+奥希替尼)把全流程走一遍、每一步都是
     真实数字——[12-worked-example-walkthrough.md](12-worked-example-walkthrough.md)
   - 想知道数据从哪来、分子怎么从骨架+片段拼出来、三路线战略的完整推导——
     [13-design-strategy-and-data.md](13-design-strategy-and-data.md)
   - 想知道一个真实药物设计师拿到"奥希替尼耐药"这个任务，每一步具体在
     想什么、下载的数据到底回答了什么问题、骨架怎么参考已有药物设计出来——
     [17-designer-thought-process.md](17-designer-thought-process.md)
   - 想看工业级差距的完整分维度表格、retrospective验证的详细结果、
     三路线各自的缺口排序——[11](11-industrial-gap-and-roadmap.md)/
     [14](14-retrospective-validation-results.md)/
     [15](15-route-specific-gaps-and-next-steps.md)
   - 想看某条路线的生物学背景——`doc/routes/*.md`
   - 想看某个环节的详细字段/阈值设计——`doc/01`~`doc/09`

## 2026-10 更新：从单一路线扩展成三条并行路线

本仓库最初只做一件事：四代EGFR TKI(克服C797S耐药)。2026-10做立项前的竞品
核查后(完整记录见[doc/competitive-landscape.md](competitive-landscape.md))，
发现这是一个竞争激烈、护城河已被部分攻破的赛道——于是仓库扩展成三条**并行**
研究的路线，不是三选一：

| 路线 | 名称 | 核心思路 | 配置目录 | 详细文档 |
|---|---|---|---|---|
| A | SHP2/SOS1变构抑制剂 | 联用EGFR TKI，打下游收敛点(MET/HER2/PIK3CA/RAS通路重新激活，合计25-30%) | `routes/route_a_shp2_sos1/` | [route-a-shp2-sos1.md](routes/route-a-shp2-sos1.md) |
| B | 突变选择性EGFR降解剂 | 催化式作用机制，同时覆盖C797S和其他靶上突变 | `routes/route_b_degrader/` | [route-b-degrader.md](routes/route-b-degrader.md) |
| C | 四代EGFR TKI(C797S) | 原始路线，最成熟但竞争最激烈 | `routes/route_c_4th_gen_tki/` | [route-c-4th-gen-tki.md](routes/route-c-4th-gen-tki.md) |

三条路线**共用同一套计算引擎**(`core/`下12个模块)，区别只在于各自
`routes/<route_id>/config/`下的配置文件(TPP表、基因型/靶点、骨架、口袋片段库)。
这个架构选择不是预先设计好的，是做完路线C之后再扩展到A/B时发现的事实：
环节1(标准化)/3(枚举)/6(ADMET)/7(合成可及性)/8(决策)/9(DMTA回传)这六个环节
完全不依赖具体靶点是什么，只有环节2(结构系综)/4(打分漏斗里和具体口袋相关的
部分)/5(选择性)的**配置**(不是代码)是路线专属的。路线B额外新增了
`core/ternary_complex.py`和`core/enumerate.py::enumerate_bifunctional()`——
因为降解剂的"三元复合物"问题和"warhead-linker-E3"三组分拼接拓扑是路线A/C
完全没有的新计算问题，不能靠配置切换解决。

## 项目背景

奥希替尼（osimertinib，三代EGFR抑制剂）获得性耐药机制是碎片化的：C797S突变
只占6-7%，MET/HER2扩增各占约10-15%，还有约50%查不出明确的基因组驱动因素。
这个碎片化格局直接决定了立项策略——"收敛"比"特异"更重要，这也是为什么
本仓库同时推进三条路线而不是押注一个最大的单一分片。完整的竞品数据和立项
推理过程见[doc/competitive-landscape.md](competitive-landscape.md)。

## 管线结构：10个环节（三条路线共用这套结构）

原设计文档(路线C最初版本，见
[四代 EGFR (C797S) 抑制剂工业级研发管线设计.md](<四代 EGFR (C797S) 抑制剂工业级研发管线设计.md>))
把管线分成环节0-9，这套结构现在是三条路线共用的框架：

| 环节 | 名称 | 代码位置 | 状态 | 详细文档 |
|---|---|---|---|---|
| 0 | 目标画像(TPP) | `routes/<route>/config/target_profile.yaml` | 配置，真实，三路线各自一份 | 本文件 + 各路线文档 |
| 1 | 分子标准化 | `core/standardize.py` | **真实可跑**，路线通用 | [01](01-standardization.md) |
| 2 | 结构系综准备 | `core/structures.py` | 部分真实+占位，路线通用 | [02](02-structure-ensemble.md) |
| 3 | 类似物枚举 | `core/enumerate.py` | **真实可跑**，含两种拓扑(固定骨架+R基团 / warhead-linker-E3三组分) | [03](03-enumeration.md) |
| 4 | 多级打分漏斗 | `core/docking.py` `core/mmgbsa.py` `core/md_stability.py` `core/fep.py` `core/covalent.py` `core/ternary_complex.py`(路线B专属) | 部分真实+占位 | [04](04-funnel-docking.md) |
| 5 | 选择性引擎 | `core/selectivity.py` | 真实数学+占位计算 | [05](05-selectivity.md) |
| 6 | ADMET/结构警示 | `core/admet.py` | **真实可跑**(有精度边界声明)，规则库全路线共享 | [06](06-admet.md) |
| 7 | 合成可及性 | `core/synthesis.py` | 部分真实(SA score)+占位 | [07](07-synthesis.md) |
| 8 | 多目标决策 | `core/decide.py` | **真实可跑** | [08](08-decision.md) |
| 9 | DMTA回传与规则治理 | `core/feedback.py` | **真实可跑**(schema+分析逻辑) | [09](09-feedback.md) |

doc/01-09这9篇文档描述的是**共享引擎**，不是某条路线专属的——里面的例子
(奥希替尼、demo骨架)来自路线C，但函数签名和逻辑对三条路线都适用。
路线专属的内容在各自的`doc/routes/*.md`里。

## 核心工程原则（贯穿所有模块，不是口号）

1. **不编造数字**。凡是需要GPU集群、MD引擎、训练好的ML模型、真实专利数据库
   这类本环境不具备的基础设施才能产出的结果，对应函数一律返回
   `{"ok": False, "reason": "..."}`，把缺什么说清楚，不返回一个"看起来合理"
   的假数值。这个原则的由来：姊妹项目 `cancer/egfr-pipline` 的受体PDBQT
   准备函数实际上是 `cp receptor.pdb receptor.pdbqt`（只是换后缀名，完全
   没有做AutoDock需要的原子类型/电荷计算），导致那个项目产出的全部对接
   分数都是 0.0 kcal/mol——一个看起来能跑、实际全错的典型反面教材。
   `core/docking.py` 的模块docstring里专门记录了这个教训。这条原则在三路线
   扩展里延续得很彻底：路线A的SHP2结构、路线B的E3配体都经过RDKit交叉核对
   (分子式/分子量匹配已知公开数据)，没有凭记忆编造化学结构；
   路线B的`core/ternary_complex.py`对无法计算的三元复合物3D预测诚实占位，
   不是硬凑一个cooperativity数值。

2. **成药性是乘法不是加法**。`core/decide.py` 的多目标打分用几何平均
   (desirability score)而不是加权和，任何一项指标为0会让总分归零，
   不允许靠某一项极端值掩盖其他项的致命缺陷。

3. **规则要有豁免清单，且豁免清单要有验证过的正负对照**。
   `config/structural_alerts.yaml`(三路线共享，不在`routes/`下面)里每条
   SMARTS规则都用真实分子验证过（奥希替尼、吉非替尼作为正常分子不应该被
   误报），验证脚本见`tests/test_structural_alerts.py`。

4. **硬门槛要慎用，弱信号不升级成拒绝规则，除非precision经过验证**。
   `core/feedback.py` 的 `rule_precision_tracking()` 实现了这条治理机制：
   一条规则要从"警告"升级成"拒绝"，必须precision>0.8且在≥3个独立化合物上
   验证过。

5. **字段schema要和设计文档的表格逐字对齐**，下游代码按字段名读取，
   不允许每个模块自己发明一套命名。`db/schema_mysql.sql`(生产环境)/
   `db/schema.sql`(本地校验用) 的每张表、`core/*.py` 里的每个dataclass/
   pydantic模型都逐字照抄了原设计文档对应小节的字段表。

6. **（新增）从第一天起就按联用设计，不要按单药设计**。一线奥希替尼耐药后
   已经不存在单药能赢的局面(amivantamab+化疗是现在的标准)。三条路线的
   `target_profile.yaml`都新增了CYP药物相互作用风险、与常用联用药物不重叠
   的毒性谱、给药方案兼容性这三项TPP指标。

## 环境与依赖

见 [`requirements.txt`](../requirements.txt)。核心依赖（RDKit/pandas/numpy/
pydantic/scipy/PyYAML/pytest）在本仓库 `.venv`（Python 3.9.6）下已验证可用。

**2026-10更新：`meeko`/`vina`/`PyRosetta`/`openmm`/`pdbfixer`/`aizynthfinder`
都已经在第二个环境`.venv310`(Python 3.10)里真实装好并跑通过**(见
`scripts/setup_docking_env.sh`，这几个都是纯pip装，不需要像vina/pyrosetta
那样处理Boost/PATH的坑)——不是"未安装"，是"主.venv装不了，需要另起一个
环境"，已经用这套环境对真实下载的PDB结构(5EHR/6LUD/6SCM/4TZ4/DUD-E的EGFR
受体)跑出了真实对接分数、真实突变建模、真实蛋白短程MD平衡(见
[02](02-structure-ensemble.md))、真实retrospective验证
(见[14](14-retrospective-validation-results.md))、真实逆合成路线搜索
(见[07](07-synthesis.md)，真实跑出奥希替尼的4步合成路线，4个起始原料全部
在ZINC真实库存里)。

**配体结合复合物的真实短程MD另外装了一个conda环境**(`openff-toolkit`需要
Python≥3.11或conda解决一套复杂依赖，pip在`.venv310`里装不了)：
`mamba install -c conda-forge ambertools openff-toolkit openmm
openmmforcefields pdbfixer`(强烈建议用mamba不要用conda classic，后者解这组
依赖可能卡在"Solving environment"几十分钟)，详见`scripts/setup_docking_env.sh`
第8步。真实跑过一次：奥希替尼(真实6LUD晶体坐标) + C797S三重突变受体，10ps
真实轨迹，配体RMSD没有发散(见[04](04-funnel-docking.md))。

仍然没有的是`openfe`(FEP，需要GPU才有实际意义，CPU上一次FEP微扰要几天到
几周)/完整PRosettaC(真正的三元复合物预测流程，需要PatchDock/完整Rosetta
C++套件这类额外的重型软件依赖)，原因和升级路径见对应模块docstring。

## 快速开始

```bash
# 建库(推荐MySQL；三条路线共用一个库，靠compounds.route_id列区分)
mysql -u root -p < db/schema_mysql.sql

# 跑测试(112个单测，覆盖所有真实可跑的模块，三条路线的配置都有独立校验)
.venv/bin/python -m pytest tests/ -v

# 路线C端到端工作流(默认route_id)
.venv/bin/python -m workflows.run_discovery_round

# 路线A端到端工作流
.venv/bin/python -c "from workflows.run_discovery_round import run; run(route_id='route_a_shp2_sos1', batch_size=10)"

# 路线B端到端工作流(枚举+ADMET画像，决策阶段诚实缺位，见该工作流docstring)
.venv/bin/python -m workflows.run_bifunctional_round

# 单独跑某个模块的演示
.venv/bin/python -m core.standardize
.venv/bin/python -m core.admet
.venv/bin/python -m core.decide
.venv/bin/python -m core.ternary_complex
.venv/bin/python -m core.route_config
```

## 已知的、故意保留的局限

- `core/admet.py` 的 `clogd_approx` 直接复用 `clogp`(中性分子分配系数)，
  不是真正按pH7.4电离态算的cLogD——需要pKa预测器才能做对，模块顶部有
  完整声明。这个局限对路线B尤其重要，因为`compute_cns_mpo()`是按传统
  小分子(MW<500)校准的，双功能分子天然beyond Rule-of-5，算出来的分数会
  系统性偏低，不代表分子"不行"，只是这把尺子不是为这类分子设计的。
- 三条路线的骨架/片段库都明确标注`illustrative_only`，不是真实候选起始点，
  仅用来演示数据结构怎么驱动`core/enumerate.py`。
- `patent_landscape.yaml`(每条路线各一份)是占位数据
  (`status: placeholder_not_real`)，不能用于真实FTO决策。
- 环节2/4(L1往后)/5/9的真实计算都需要本环境不具备的基础设施(PDB结构库+
  MD软件+GPU集群、湿实验数据)，当前是正确接口+诚实占位。路线B额外缺失
  三元复合物3D预测工具(PRosettaC/Rosetta)，接入方式见各自的分文档。
