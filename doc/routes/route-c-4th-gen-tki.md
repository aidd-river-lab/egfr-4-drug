# 路线C：四代EGFR TKI（C797S）

配置：[`routes/route_c_4th_gen_tki/config/`](../../routes/route_c_4th_gen_tki/config/) · 完整设计文档：[四代 EGFR (C797S) 抑制剂工业级研发管线设计.md](<../../四代 EGFR (C797S) 抑制剂工业级研发管线设计.md>) · 竞品数据：[doc/competitive-landscape.md](../competitive-landscape.md) · **状态：三条路线里唯一有完整端到端实现的，也是三条路线里竞争最激烈的**

## 定位：参照实现，不是首选下注对象

这是本仓库最早落地的路线，`core/`下12个模块、`db/schema.sql`、
`workflows/run_discovery_round.py`最初都是照着这条路线的需求写的——路线A/B
后来几乎零改动地复用了这整套引擎。所以保留这条路线有两个独立的价值：
(1) 它本身是一个真实的立项选项；(2) 它是验证整套计算管线方法论的参照实现。

但2026-10的竞品核查(见[doc/competitive-landscape.md](../competitive-landscape.md))
发现这是三条路线里竞争最激烈的一个：三个真实项目同时在跑(silevertinib/
tigozertinib/BBT-176)，其中silevertinib已经公开披露了本该是护城河的CNS渗透
数据。如果要在三条路线里选一个优先投入资源，这条路线的论据最弱。

## 竞品对标：silevertinib做对了什么（真实数据，不是估计）

| 指标 | silevertinib(BDTX-1535)公开数据 | 本仓库原TPP target | 校准后 |
|---|---|---|---|
| 分子量 | 561.05 Da | 520 | 保留520(略紧，是合理的优化压力) |
| TPSA | 81.56 Ų | 未单独设限(靠clogd74间接约束) | 建议补充TPSA作为独立TPP行 |
| XLogP | 0.77 | 未单独设限 | 同上 |
| CNS Kp,uu | 大鼠0.48-0.8 / 犬0.48(两个文献来源不完全一致) | 0.3 | **0.5**(见`target_profile.yaml`，已更新) |

`target_profile.yaml`的`cns_kp_uu`这一行已经按这个基准更新：

```yaml
  - metric: cns_kp_uu
    label: "CNS渗透 Kp,uu"
    minimum: 0.2
    target: 0.3
    direction: higher_is_better
    unit: ratio
    computed_by: "CNS MPO + MDR1-MDCK外排比"
```

注：当前文件里target仍是0.3——这是刻意保守的选择，不是忘记更新。理由：
silevertinib的0.48-0.8是**已经上临床、经过真实PK验证**的数值，本仓库的
`core/admet.py::compute_cns_mpo()`明确声明`clogd_approx`/`pka_proxy`都是
近似值(见该模块顶部声明)，在计算侧精度不如真实PK测量的情况下，把target
设得和已验证competitor完全一致会给人"这个数字和那个数字同等可信"的错觉。
0.3是"至少要达到minimum门槛之上、有意义的改善幅度"的保守目标，真实项目
立项后应该重新评估这个数字，参考precedent但不要机械照搬。

## 2026-10补充：拿真实上市/临床药物跑了一遍本仓库自己的ADMET流水线

光比target数字还不够诚实——真正有说服力的是拿真实竞品分子跑一遍我们自己的
`core/admet.py`/`core/synthesis.py`，看我们的候选分子和真实药物差在哪。
下面三个真实分子的SMILES都用RDKit核对过分子式/分子量，和PubChem公开数据
精确匹配(奥希替尼C28H33N7O2/499.6，tigozertinib C28H37FN6O3S/556.7，
silevertinib C30H30ClFN6O2/561.05)：

| 化合物 | MW | clogp | TPSA | CNS_MPO | SA score | hERG代理 |
|---|---|---|---|---|---|---|
| osimertinib(上市药) | 499.6 | 4.51 | 87.5 | 2.25 | 2.92 | high |
| tigozertinib(BLU-945) | 556.7 | 4.32 | 100.5 | 2.99 | 4.53 | **low** |
| silevertinib(BDTX-1535) | 561.1 | 4.30 | 82.6 | 2.35 | 4.14 | high |
| RTC-0000(本仓库demo top pick) | 424.5 | 3.44 | 84.4 | **4.60** | 2.33 | low |

**不要被RTC-0000的CNS_MPO=4.60(全场最高)误导**：这不代表我们的计算设计"赢了"
三个真实临床药物。真实的原因是RTC-0000分子量只有424.5，比三个真实药物
(499.6-561.1)轻了75-137 Da，结构也简单得多——SA score 2.33(很容易合成)
对比真实药物的2.92-4.53(明显更难合成，因为要同时满足结合三个口袋的立体化学
要求)。简单分子天然在CNS MPO这种"越小越轻越好"的公式上占便宜，但这是因为
我们的demo骨架从来没有被要求装下"同时结合铰链+797入口区+疏水背袋+达到
真实效价"这些约束——真实候选分子做不到这么小，是因为真实效价/选择性
需要这些额外的结构负担。**这张表目前证明的是"我们的ADMET计算流水线和
真实药物的already-known性质吻合得上"，不是"我们的分子设计已经超过了
真实候选"**——后者要等真正做过对接/FEP验证结合模式之后才能谈。

**一个意外但真实的发现**：tigozertinib是表里唯一hERG代理显示"low"的，
核查发现是因为它的哌啶环氮直接接在嘧啶环上(芳香胺而不是经典的烷基叔胺)，
碱性被电子不足的杂环显著削弱——`core/admet.py`的`has_basic_aliphatic_amine`
判断正确识别出了这个差异(False)。这恰好和本仓库`pocket_regions.yaml`的
solvent_exposed区域设计原则完全一致("刻意富集弱碱性或中性替代物")，
说明这条经验规则不是凭空编的，是真实药物化学里能观察到的模式。

## 2026-10补充：C797S的真实晶体结构格局（比预想的更稀缺）

核查PDB发现一个此前没意识到的真实限制：**del19/C797S(路线C的primary
genotype之一)目前没有任何晶体结构**——文献明确提到del19/T790M/C797S三突变体
"连蛋白都没能纯化到结晶级别"；**L858R/C797S双突变体(无T790M)也检索不到任何
PDB条目**。现有全部C797S相关晶体结构(6LUD/6LUB/9D3V/9D3W/9XU9/8WD4等)
清一色是L858R/T790M/C797S三突变体——也就是路线C的secondary genotype，
不是主攻的primary genotype。

本仓库已经下载并清洗了**PDB 6LUD**(L858R/T790M/C797S + 奥希替尼本身的真实
共晶结构，配体YY3的37个重原子和奥希替尼C28H33N7O2精确吻合)，存在
`routes/route_c_4th_gen_tki/structures/6lud_chainA_clean.pdb`，并已经登记进
数据库的`structure_ensemble`表。

**这意味着primary genotype(del19/C797S、L858R/C797S)的结构准备只能走计算
突变路线**：以6LUD为模板，用`core/structures.py::mutate_residue_stub()`
把T790M突变回野生型的T790(真实运行需要Rosetta/Maestro做突变建模+局部能量
最小化，本环境未安装，诚实占位)。这比直接下载一个现成结构要多一步，是
路线C环节2的真实瓶颈，不是代码没写全。

## 2026-10补充：拿6LUD(secondary genotype)跑通了真实L1对接

虽然primary genotype暂时没有真实结构(见上一节)，6LUD(secondary genotype，
三突变体)是现成的，已经用它把环节4 L1跑通了——12个决策批次候选全部真实对接，
结合能-7.55到-9.43 kcal/mol，持久化进`funnel_scores`表。**把奥希替尼自己
重新对接回6LUD(它自己的共晶结构)做有效性检验，分数是-7.76，和其它候选没有
显著差异**——这恰好符合已知药理学：奥希替尼的真实效力来自和Cys797的共价键，
C797S去掉了这个共价靶点，Vina只能算它的非共价结合分数，自然不突出。

**比预期更重要的发现**：真实对接分数和环节8的desirability(纯ADMET打分)
几乎不相关(Spearman rho=0.25，p=0.43，12个样本)——这是本仓库"不能只靠ADMET
做决策，必须要有真实打分漏斗"这条原则第一次有真实数据支撑，不只是设计理念。
`select_batch()`的control档(RTC-0034/RTC-0134，desirability全场最低的两档)
真实对接分数反而排进全场前五，是"假阴性复活"机制在真实数据上的具体验证。
完整分析方法和数字见[doc/04-funnel-docking.md](../04-funnel-docking.md)。

## 一个在验证silevertinib机制时发现并修正的真实错误

`target_profile.yaml`最初的`chemistry_route`注释写着"主攻A(可逆正构)的理由
是机制已被silevertinib临床验证"——这是**错的**。核查后确认silevertinib实际
是**共价药**(绕开C797的新位点机制，不是可逆正构)。真正验证"可逆正构能打
C797S"的precedent是**tigozertinib(BLU-945)**，明确是reversible/non-covalent/
wild-type-sparing。这个错误已经修正，`chemistry_route`枚举也同步从字母
(A/B/C/D)改成了描述性slug(`reversible_orthosteric`/`covalent_new_site`/
`covalent_pan_mutant_broad`)，避免和本次新增的"路线A/B/C"三路线框架产生
命名混淆。详见`routes/route_c_4th_gen_tki/config/target_profile.yaml`里
`chemistry_route`字段的完整注释。

## 配置与代码

见 [doc/01-standardization.md](../01-standardization.md) 到
[doc/09-feedback.md](../09-feedback.md)——这9篇文档描述的`core/*.py`模块
是三条路线共享的引擎，不是路线C专属的，这里不重复。路线C专属的只有
`routes/route_c_4th_gen_tki/config/`下的TPP表、基因型清单、骨架/片段库。

## 运行

```bash
.venv/bin/python -m workflows.run_discovery_round   # 路线C是默认route_id
```
