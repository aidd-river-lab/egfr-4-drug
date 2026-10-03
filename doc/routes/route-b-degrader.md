# 路线B：突变选择性EGFR降解剂

配置：[`routes/route_b_degrader/config/`](../../routes/route_b_degrader/config/) · 新增代码：[`core/ternary_complex.py`](../../core/ternary_complex.py) · 竞品数据：[doc/competitive-landscape.md](../competitive-landscape.md) · **状态：枚举+ADMET画像可跑，决策阶段诚实缺位（见下）**

## 为什么是这条路线

和路线C覆盖同样的病人群(C797S等靶上耐药突变)，但作用机制完全不同：不是
占据ATP口袋抑制激酶活性，而是招募E3连接酶把突变EGFR蛋白整个降解掉。
四个理论优势：
1. **催化式作用机制**，低暴露量即可起效，绕开四代TKI最难的"要压过细胞内
   毫摩尔级ATP浓度"这个问题；
2. 理论上能**同时覆盖C797S和其他靶上突变**——降解不挑具体突变类型，
   只要warhead还能结合该突变体就行；
3. 顺带去掉EGFR的激酶非依赖支架功能；
4. 耐药需要E3连接酶本身丢失，是更罕见的逃逸路径。

真实precedent是CFT8919(C4 Therapeutics，EGFR L858R选择性降解剂)——核查后
发现它比想象中更早期：IND刚获批，中国1期试验"尚未开始入组"(见
[doc/competitive-landscape.md](../competitive-landscape.md))。这意味着整个
赛道几乎没有临床验证数据：技术壁垒高(护城河)，但风险也真实存在。

## 最大的障碍：CNS渗透

这不是本仓库编的风险——文献明确指出"大多数EGFR降解剂体外活性很好，但PK
性质(尤其CNS穿透)往往跟不上"。所以`target_profile.yaml`把`cns_kp_uu`
提升为一级TPP指标(不是像路线C那样的常规监测项)，target值(0.3)故意设得
比路线C的校准值(约0.5-0.8，对标silevertinib)更保守——因为双功能分子天然
beyond Rule-of-5，这个目标本身已经是本路线最难达成的指标。

## 技术新增：三元复合物问题

二元对接(warhead结合EGFR、E3配体结合CRBN)只是必要条件，不是充分条件——
真正决定降解效率的是**三元复合物**(warhead-EGFR-CRBN)能不能形成，以及
形成时两边的结合界面是互相促进(正协同)还是互相干扰(负协同)。这是二元
对接完全捕捉不到的，需要PRosettaC/Rosetta这类专门工具，本环境未安装。

`core/ternary_complex.py`延续"真实数学+诚实占位"的模式：
- `run_ternary_complex_stub()`：诚实返回`ok=False`，不编造cooperativity数值。
- `evaluate_hook_effect_risk()`：**真实、可测试的部分**。Hook effect(钩状效应)
  是双功能降解剂公认的真实药理学现象——高浓度下warhead和E3配体各自的二元
  结合都趋于饱和，反而会让三元复合物比例下降，剂量-响应曲线呈钟形。这条
  判据只需要两个二元Kd和给药浓度范围，不需要三元复合物3D预测就能算，
  已用6个场景测试覆盖(`tests/test_ternary_complex.py`)。

## 新增的枚举拓扑：warhead-linker-E3三组分

和路线A/C"固定骨架+R基团"的拓扑不同，双功能分子没有"谁是骨架"的区分——
三个独立片段依次相连。`core/enumerate.py::enumerate_bifunctional()`复用
同一套molzip+连通性校验逻辑(这类bug和骨架拼接是同一类风险，甚至更容易出错，
因为是两个独立的同位素配对而不是单侧配对)。

```yaml
warhead_library.yaml   # 复用路线C已验证的EGFR结合片段，衍生出1个连接点的warhead
linkers.yaml            # PEG链(peg1/2/3)和三氮唑，2个连接点
e3_ligands.yaml          # 来那度胺/泊马度胺的N-连接衍生物，1个连接点
```

E3配体的真实性已核实：来那度胺(C13H13N3O3，MW 259.26)、泊马度胺
(C13H11N3O4，MW 273.24)都是已上市药物，用RDKit交叉核对过分子式与公开数据
库一致。**VHL配体(如VH032)没有收录**——VH032含hydroxyproline骨架和多个
立体中心，本仓库没有核实过精确SMILES，为避免凭记忆复现出立体化学错误的
结构，诚实地留空而不是编一个"看起来像"的近似结构。

1个warhead × 5个linker × 2个E3配体 = 10个候选，全部验证连通
(`tests/test_route_b_rgroup_libraries.py` + `tests/test_enumerate_bifunctional.py`)。

## 运行

```bash
.venv/bin/python -m workflows.run_bifunctional_round
```

## 这个工作流诚实缺失的部分

和路线A/C的`run_discovery_round.py`不同，`run_bifunctional_round.py`**没有
跑多目标决策(环节8)**——因为路线B的核心TPP指标(dc50_nm/dmax_pct/
hook_effect_margin)都需要warhead/E3配体各自的结合常数(Kd)作为输入，而这些
Kd要么来自真实FEP/对接(本仓库诚实占位)，要么来自已发表的E3配体文献值
(可以查到，但本次没有接入)。**没有真实Kd就不编一个假的去跑决策算法**——
接入路径见`workflows/run_bifunctional_round.py`文件末尾的TODO注释。
