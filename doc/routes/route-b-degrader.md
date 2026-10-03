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

## 2026-10补充：E3配体的真实Ki找到了，warhead的Kd还没有

查文献确认了来那度胺/泊马度胺对CRBN的真实结合常数：**来那度胺Ki≈3.1 μM
(3100 nM)，泊马度胺Ki≈0.8 μM(800 nM)**(不同assay格式下数值有一定差异，
这在该领域是常态，不是数据矛盾)。这意味着`core/ternary_complex.py::
evaluate_hook_effect_risk()`现在可以用**真实的E3配体侧Kd**跑：

```python
evaluate_hook_effect_risk(warhead_kd_nm=100, e3_ligand_kd_nm=3100, dose_range_nm=(10, 1000))
```

但`warhead_kd_nm`这一侧依然没有真实值——它需要warhead对突变EGFR的真实结合
数据，来自环节4的对接/FEP(本环境工具链未装，诚实占位)，查文献也查不到
(因为我们的warhead是从路线C demo骨架衍生出来的illustrative结构，不是
已发表化合物)。**这是目前唯一真正卡住路线B决策阶段的缺口**：E3配体那一半
数据是真的，warhead那一半还是假设值，hook effect可以演示逻辑但不能用于
真实决策。

## 2026-10补充：三元复合物几何筛选真实跑了一次，找到一个真实的负向发现

放弃完整复刻PRosettaC发表论文的确切流程(需要PatchDock+完整Rosetta C++套件+
PyMOL+HPC调度系统，这些在单台Mac上不现实)，改用PyRosetta自带的刚体对接协议
(`RigidBodyRandomizeMover` + `FaDockingSlideIntoContact`，不需要PatchDock/
调度系统)搭了一个更轻量的几何筛选版本：

1. 真实蛋白-蛋白对接：EGFR(6LUD来源) vs CRBN(4TZ4来源)，两条链都是标准氨基酸，
   不需要把warhead/来那度胺参数化成Rosetta自定义残基类型(这条路线试过，
   碰到`RDMolToRestype`需要先转换成PyRosetta自己编译的RDKit绑定对象，
   比标准Python rdkit.Chem.Mol多一层转换，文档稀少；后来改用Rosetta原生的
   `SDFParser`+`convert_to_ResidueType`路径验证能解析分子，但构建完整可用的
   自定义ligand residue还需要更多工作——这条路线最终放弃，不是做不到，是
   工程代价远超当前阶段需要)。
2. 用`core/ternary_complex.py::kabsch_rigid_transform()`(纯numpy的Kabsch算法)
   在"蛋白链对接前后CA坐标的变化"里反推出刚体变换，再把这个变换套用在
   配体的exit vector原子上——这样完全不需要把配体做成Rosetta能懂的残基，
   只用蛋白骨架的坐标就能追踪配体跟着挪到哪了。
3. warhead的exit vector：真实对接进EGFR口袋拿到的3D坐标(用碘原子临时标记
   连接点，docking后按元素唯一定位，不会被原子重排搞乱)。
4. 来那度胺的exit vector：4TZ4真实晶体结构里的N17原子(靠原子间距离推断连接关系
   确认的，不是读PDB原子名瞎猜的——N17只连了一个重原子C14，是经典的末端
   芳香胺模式，和来那度胺的4-氨基异二氢吲哚酮结构精确吻合)。
5. 每根linker的真实span：用RDKit生成50个MMFF优化构象，测连接点间距离分布
   (不是凭经验猜的数字)——peg1 2.5-5.9Å，peg2 4.3-8.7Å，peg3 6.7-11.9Å，
   triazole_short 5.6-8.1Å，alkyl_c4 2.0-4.6Å。

**真实跑出来的结果**：500次独立随机刚体对接+滑入接触，exit vector距离最小
只有26.8Å，中位数56.4Å——**0/500次落进任何一根linker的真实span范围内**。

这是一个真实、有统计量支撑的负向发现，不是bug：均匀随机的刚体朝向采样，
两个蛋白表面上各自很窄的exit vector朝向恰好对上的概率天然很低。这正是
真实PRosettaC论文要用PatchDock(系统性穷举表面形状互补的补丁，不是均匀随机
采样)而不是随机对接的原因——本仓库这次没有完整复刻那条路径，验证到的是
"轻量版行不通"这个真实的工程结论，而不是"三元复合物不存在"。

要让这条路线真正跑通，下一步大概率需要：
- 用已知的exit vector大致方向(而不是整个蛋白表面)去约束初始刚体朝向的采样范围，
  不再是均匀随机
- 或者接受更重的方案(真的装PatchDock——需要解决macOS兼容性；或者用已经装好的
  PyMOL+现有工具自己写一个表面形状互补的候选姿态生成器)

`core/ternary_complex.py`里新增的`run_protein_protein_docking_trials()`/
`evaluate_ternary_complex_geometry()`/`measure_linker_span()`/
`kabsch_rigid_transform()`都是可复用的真实函数(lazy import PyRosetta，
纯数学部分不需要)，不是这一次性脚本——下次换更好的采样策略时，这些函数
都能直接复用。

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
