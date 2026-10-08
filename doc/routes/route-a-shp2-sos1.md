# 路线A：SHP2/SOS1变构抑制剂（联用EGFR TKI）

配置：[`routes/route_a_shp2_sos1/config/`](../../routes/route_a_shp2_sos1/config/) · 竞品数据：[doc/competitive-landscape.md](../competitive-landscape.md) · **状态：配置完整，和路线C共用同一套core/引擎**

## 为什么是这条路线

奥希替尼耐药后约25-30%的病人是通过下游RAS-MAPK通路重新激活逃逸的
(MET扩增、HER2扩增、PIK3CA突变、其他旁路重新激活)——这些机制各自独立看
都是小分片，但它们最终都收敛到同一个下游节点。打这个收敛点，理论上能用
一个联用搭档同时覆盖好几个分片，而不是每个机制单独开发一个药。

这和路线C(四代EGFR TKI)的立项逻辑完全不同：路线C是"更好地打EGFR本身"，
路线A是"不管EGFR耐药的具体分子机制是什么，在下游拦一道"。

## 为什么这条路线对计算化学友好

SHP2的变构tunnel位点结构表征非常充分——PDB **5EHR**(SHP099复合物)是这个领域
的经典参考结构，关键相互作用残基(Phe113/Arg111/Glu250)清楚，加上TNO155/
RMC-4630等临床化合物的公开数据，FEP相对结合自由能在这个口袋上能发挥作用
(`core/fep.py`)。这是`routes/route_a_shp2_sos1/config/target_profile.yaml`
把路线A排在"计算化学能发力程度：高"的原因。

但要诚实地说：**这不等于"容易"**。核查发现TNO155单药在118例重度预处理患者中
仅20%达到疾病稳定——靶点生物学验证了，但"联用策略本身是否有效"这件事还没
被充分证实。真正的挑战在SHP1选择性(最近的同源旁系基因，tunnel位点高度相似)
和联用场景下的真实疗效，不在"能不能算出好的对接分数"。

## 配置结构

- `target_profile.yaml`：SHP2(tunnel位点)为主攻，SOS1(催化邻近口袋)为backup；
  可及人群按**耐药机制分型**(不是EGFR突变分型)；TPP里DC50/Dmax不适用(这不是
  降解剂)，核心终点是SHP2生化IC50+SHP1选择性倍数+联用细胞活性。
- `pipeline.yaml`：和路线C一样的L0-L4打分漏斗形状，`hinge_hbond`/`target_anchor`
  这两个key被挪用表示tunnel位点的关键H-bond网络占据率(不是激酶铰链/Ser797)——
  复用`core/md_stability.py::evaluate_pass_criteria()`不需要改代码，只是语义
  在pipeline.yaml的注释里重新定义了。
- `rgroup_libraries/`：`demo_aminopyrazine_tunnel`骨架，灵感来自已发表的
  "aminopyrazines as allosteric SHP2 inhibitors"化学型文献，**不是对SHP099或
  任何specific化合物的精确重现**(避免凭记忆复现复杂结构出错)。4个tunnel胺臂
  片段 × 4个远端联芳基片段 = 16个组合，已验证全部连通(见
  `tests/test_route_a_rgroup_libraries.py`)。

## 运行

```bash
.venv/bin/python -c "
from workflows.run_discovery_round import run
run(route_id='route_a_shp2_sos1', batch_size=10)
"
```

跑的是和路线C同一个工作流(`workflows/run_discovery_round.py`)，只是
`route_id`参数不同——这是"共享引擎+路线专属配置"架构的直接体现：枚举/
标准化/ADMET/决策四个环节一行代码不用改。

## 和路线C的关键差异

| | 路线A | 路线C |
|---|---|---|
| 可及人群分型依据 | 耐药机制(MET扩增等) | EGFR突变(C797S等) |
| 化学机制 | 非共价变构 | 可逆正构/共价(待定) |
| 主终点 | 生化IC50 + 选择性倍数 | 细胞IC50 + WT选择性 |
| 必须联用 | 是，从立项起就是联用设计 | 可以单药，但联用也是大趋势 |
| 复用的core/模块 | 100%复用，无需新增代码 | 原始实现 |

## 2026-10补充：真实结构已下载

`5EHR`(SHP099/SHP2，tunnel位点)和`6SCM`(BI-3406/SOS1，催化邻近口袋)已经
实际下载+清洗(保留了共晶配体)，存在`routes/route_a_shp2_sos1/structures/`，
并登记进数据库的`structure_ensemble`表。确认了配体HETATM确实存在
(5EHR里的5OD=SHP099，23个重原子；6SCM里的L7H=BI-3406，**33个重原子**——
之前这里写错成"58个重原子"，核实后发现58是把PDB里25个显式氢原子也算进去
的总原子数，真实重原子数是33，已修正)，
不是空结构。路线A这一步比路线C更省力：SHP2/SOS1的真实共晶结构本身就
存在，不需要像路线C的primary genotype那样走计算突变建模(见
`doc/routes/route-c-4th-gen-tki.md`的说明，那边连晶体结构都没有，只能靠
PyRosetta计算产出一个近似模型)。

## 2026-10补充：L1真实对接已跑通，16个候选全部有真实结合能

用`scripts/setup_docking_env.sh`搭的`.venv310`对全部16个demo候选跑了真实Vina
对接(受体=5EHR，口袋中心取自SHP099的真实坐标质心)，结合能范围
-7.73到-10.74 kcal/mol，全部持久化进`funnel_scores`表(funnel_level=L1)。
**发现**：这批对接分数和环节8算出来的desirability(纯ADMET打分)的Spearman
相关系数是0.83(p=0.04，6个样本)——比路线C的同类分析(见
[doc/04-funnel-docking.md](../04-funnel-docking.md))相关性更强，但样本量太小
(n=6)不能太当真。更可靠的信号是：`select_batch()`的control档(RTA-0005，
desirability全场最低)真实对接分数反而是全场最佳(-10.74)——"假阴性复活"
机制在真实数据上又验证了一次。

## 2026-10再补充：6SCM(SOS1)真实准备好了，但16个候选从来没对接过去

之前只对5EHR(SHP2)跑过真实对接，6SCM(SOS1)的受体PDBQT从来没准备过——
现在补上了：真实清洗+meeko准备(`--allow_bad_res --default_altloc A`，
6SCM里9个残基有alternate location，和5EHR/6LUD处理方式一致)。做了一次
**自我一致性验证**(和之前奥希替尼redocking进6LUD/4ZAU同一个思路)：
把BI-3406自己真实的SMILES(查证过，ChEMBL CHEMBL4519023，formula
C23H25F3N4O3/MW 462.5，匹配真实发表数据)重新对接回6SCM，box中心取自
真实晶体坐标质心——**真实打分-8.95 kcal/mol，姿态空间延展10.2×11.1×8.4Å**
(合理的单分子尺度，不是之前路线B那次"整分子被压扁"的那种假象)。这个
分数量级和BI-3406真实发表的IC50=6nM(强效)方向一致，是一次健康的sanity
check。

**顺带纠正了一个之前记录错误的数字**：上面"2026-10补充：真实结构已下载"
那段原来写"L7H=BI-3406，58个重原子"——核实后发现58是把PDB里25个显式氢
也算进去的总原子数，真实重原子数是33(已在上面修正)。

**但16个真实候选分子依然没有对接进6SCM过**——这16个候选的骨架是照着
SHP099(SHP2 tunnel位点)的化学逻辑设计的，不是照着BI-3406(SOS1催化邻近
口袋)设计的，直接拿去对接6SCM意义有限(大概率因为化学型不对口而打分不好，
不能说明"这16个候选不行"，只能说明"用错了口袋"——和路线A那次SHP2
ChEMBL验证踩到的坑是同一类问题)。如果真要评估这16个候选对SOS1有没有
交叉结合风险(selectivity检查的一部分)，这个真实对接是有意义的下一步，
但目前没有做。

## 2026-10再补充：SOS1专属retrospective验证真实跑完了——结果比SHP2和EGFR都强

上面那次只是BI-3406单分子的sanity check，不是"这把尺子在SOS1上准不准"
的系统性验证。补了一次和SHP2验证同结构的真实ChEMBL检验：查CHEMBL2079846
(人源SOS1)真实IC50/Ki，active(≤1μM)/inactive(≥10μM)筛出776/19个
(负例比SHP2那次的279个少得多，是SOS1这个靶点比SHP2年轻、历史文献里
"真测过但无效"的记录天然更少，不是抽样问题)，active抽样400个+inactive
全部19个，共419个真实化合物，对接进这次刚准备好的6SCM受体(box中心用
的就是上面BI-3406 sanity check验证过的真实晶体质心)。

**真实结果，419个全部成功**：

| 指标 | SOS1(6SCM) | SHP2(5EHR) | EGFR(DUD-E) |
|---|---|---|---|
| AUC | **0.878** | 0.575 | 0.69 |
| active/inactive平均打分 | -9.10 / -7.66 | -8.42 / -8.41 | -9.70 / -8.89 |

AUC=0.878是这个仓库三次retrospective验证里最高的一次，分子量混淆因素
排查过(active分子量反而比inactive小，打分-分子量相关系数仅0.094)，
排除了"分数只是在追踪分子大小"这个可能性。EF指标这次不拿来当主要
依据——active占了95.5%的极端类别比例下，EF的理论上限本身就逼近1，
不是方法失效，是这个场景下EF没有信息量，完整解释见
[doc/15](../15-route-specific-gaps-and-next-steps.md#2026-10再补充路线a的sos1专属retrospective验证真实跑完了结果比shp2和egfr都强)。

**这件事改变了一个之前过于笼统的判断**：之前"路线A不如EGFR可信"是
拿SHP2那次的弱结果代表整条路线，但SOS1这次明显更强——同一套Vina
代码在路线A内部两个受体上的表现差了一大截，不能一概而论。如果要在
路线A内部选一个受体优先投入下一级验证(MD/FEP)，这次真实数据支持
优先选SOS1。注意这**不等于**"16个候选分子对SOS1真的有效"——这次验证
用的是ChEMBL真实已知活性的外部化合物，不是本路线自己设计的16个候选
(那16个还是照SHP099化学型设计的，见上一节，依然没有对接进6SCM)，
"对接方法在这个受体上区分力强"和"我们设计的候选恰好命中这个受体"是
两件独立的事，不能混为一谈。

原始数据：`validation/chembl_sos1/results/docking_scores.csv`；脚本：
`scripts/run_chembl_sos1_validation.py`。

## 2026-10再补充：真实在免费GPU上跑完了3个真实候选分子的复合物MD——整个仓库第一次

`notebooks/colab_gpu_md.ipynb`之前一直没人真实在GPU上跑过。这次用
`scripts/export_top_candidates_for_gpu_md.py`导出了路线A(SHP2/5EHR)
真实L1对接分数最好的3个候选(RTA-0005/-10.74、RTA-0009/-9.5、
RTA-0010/-9.85 kcal/mol)，真实在Colab免费T4 GPU上跑完了复合物MD
(走的是OpenCL平台，不是CUDA——过程中CUDA插件真实遇到了PTX版本不兼容，
notebook已经改成自动探测能用的平台，详见`doc/04`的debug记录)。

**真实结果**(0.474ns/候选，20帧，7865-7875原子，GPU吞吐6.6秒/1000步，
比本地CPU的17.2秒/1000步快约2.6倍)。轨迹趋势用一个朴素但可复现的量化
方法判定(前半段10帧均值 vs 后半段10帧均值，差超过0.5Å才算有明显趋势，
不是凭眼睛看形状)：

| 候选 | 均值RMSD | 前半段→后半段 | 趋势 | L1对接分数 | R1基团 |
|---|---|---|---|---|---|
| RTA-0005 | 3.33Å | 4.0Å → 3.0Å | **settling**(收敛) | -10.74(三个里最好) | 4-aminopiperidin-4-yl |
| RTA-0009 | 3.61Å | 3.19Å → 4.03Å | **drifting**(漂移) | -9.5 | 3-aminopyrrolidin-1-yl |
| RTA-0010 | 3.27Å | 2.95Å → 3.6Å | **drifting**(漂移) | -9.85 | 3-aminopyrrolidin-1-yl |

**老实地说，三个均值都没有达到`pipeline.yaml`里L3的真实pass门槛(<2.5Å)**，
但这次轨迹只有0.474ns，离真实项目要求的15ns窗口差了三个数量级，均值
本身不能直接套用pass criteria下结论——**更有信息量的是趋势**：唯一
settling(收敛)的RTA-0005，恰好是L1打分最好、且R1基团是
`4-aminopiperidin-4-yl`的那个；两个drifting(持续漂移，更像真实不稳定，
L1给的对接姿态可能是假阳性)的都用了`3-aminopyrrolidin-1-yl`这个R1。
**这是一个假设，不是结论**——样本量只有1 vs 2，远不够格说"这个R1基团
不行"，但值得在剩下13个候选里重点验证：这条路线demo骨架的R1位点一共
4个选项(`4-amino-1-methylpiperidin-4-yl`/`4-aminopiperidin-4-yl`/
`3-aminopyrrolidin-1-yl`/`3-aminoazetidin-1-yl`，各4个R2搭配)，目前
一个`4-aminopiperidin-4-yl`+两个`3-aminopyrrolidin-1-yl`已经测过，
另外两种R1(甲基化哌啶/氮杂环丁烷)一个都没测过。

**两个drifting的候选已经记录了结构化假阳性归因**(`false_positive_
attributions`表，标签`pose_error`，`attributed_by`老实标成
`automated_rmsd_trend_heuristic_v1`，不冒充人工复盘)——但`core/feedback.py`
的治理原则要求precision>0.8且≥3个独立化合物才能把假设升级成
`structural_alerts.yaml`里的硬拒绝规则，现在连"3个独立化合物"都不够，
**这次只是开始积累数据，不是已经找到了可以拒绝的R基团**。

真实数据已经落库(`funnel_scores`表L3级+`md_stability_qc`表，
`qc_pass`老实留NULL；`false_positive_attributions`表记了上面两条
pose_error归因)，脚本：`scripts/persist_gpu_md_results.py`；原始轨迹：
`notebooks/my_candidates_gpu_md_results/RTA-xxxx_md.dcd`。

### 2026-10再补充：剩下13个候选已经导出，等着在GPU上验证上面那个R1假设

`scripts/export_top_candidates_for_gpu_md.py`现在默认跳过已经真实跑过
L3的候选，重新跑一次`--route route_a_shp2_sos1 --top-n 16`真实导出了
剩下13个(`RTA-0006/0011/0008/0007/0012/0001/0004/0014/0013/0015/0002/
0003/0000`)，打包在`validation/gpu_md_export/route_a_shp2_sos1_top16.zip`，
已经带上每个候选的真实R1/R2基团组成(写进了MANIFEST.txt)，方便跑完之后
直接按R1分组对比趋势，不用再回头查数据库。这13个里`4-amino-1-methyl
piperidin-4-yl`(4个)和`3-aminoazetidin-1-yl`(4个)这两种R1目前一个
真实MD结果都没有，是验证上面假设最关键的缺口。

## 尚未做的事

- L2精细重打分/共价对接(本路线不涉及共价，但可以用Gnina CNN rescoring)
- 用TNO155/RMC-4630的已发表活性数据做FEP map的实验锚点校正
- CYP/毒性谱/给药方案这三项联用TPP指标目前只有qualitative占位，需要接入
  真实CYP抑制预测模型
- 16个候选对SOS1(6SCM)的交叉结合风险检查(受体已经准备好了，差的是
  真正跑这16个候选)
