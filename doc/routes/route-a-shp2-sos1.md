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

## 尚未做的事

- L2精细重打分/共价对接(本路线不涉及共价，但可以用Gnina CNN rescoring)
- 用TNO155/RMC-4630的已发表活性数据做FEP map的实验锚点校正
- CYP/毒性谱/给药方案这三项联用TPP指标目前只有qualitative占位，需要接入
  真实CYP抑制预测模型
- 16个候选对SOS1(6SCM)的交叉结合风险检查(受体已经准备好了，差的是
  真正跑这16个候选)
