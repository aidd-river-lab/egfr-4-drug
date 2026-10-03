# 竞争格局核查记录（三路线立项决策的事实依据）

本文档记录2026-10做三路线立项决策时核实过的全部关键数据点，每一条都标注了核查
方法和来源。**后续讨论/文档/代码注释引用这些数据时，应该引用这份文档而不是凭印象
重复**——数据会过时，但"在哪验证过、验证方法是什么"这个记录本身不会过时。

## 一、一线奥希替尼耐药后的竞争格局

| 方案 | 覆盖机制 | 关键数据(已核实) | 状态 |
|---|---|---|---|
| amivantamab + 化疗 | 机制不可知 | MARIPOSA-2：确认ORR 64%(amivantamab+化疗)/63%(三药) vs 36%(化疗)；PFS风险下降56%(三药)/52%(两药) | 已获批，当前事实标准 |
| osi + savolitinib | MET扩增/高表达 | SACHI(三期)：ORR 58% vs 34%，mPFS 8.2 vs 4.5个月(ITT)或9.8 vs 5.4个月(一二代TKI后亚组)；NDA已报FDA | 接近获批 |
| HER3-DXd(patritumab deruxtecan) | 机制不可知 | HERTHENA-Lung02：中位OS 16.0 vs 15.9个月，**HR 0.98**，BLA已撤回(尽管PFS此前显著改善) | 失败 |
| 四代EGFR TKI | C797S(6-7%) | 见下方专节 | 1/2期，三个真实项目在跑 |

**核查方法**：WebSearch检索一手医学新闻源(OncLive/TargetedOnc/CancerNetwork等)和
官方新闻稿(J&J/Daiichi Sankyo/Merck)，交叉确认关键数字。

**注意和早期讨论稿的出入**：MARIPOSA-2的ORR数字此前错误引用成"53% vs 29%"，
核查后确认应为"64%/63% vs 36%"——如果这份数据要用于正式立项材料，以这份文档
为准，不要用更早的版本。

### HER3-DXd失败的含义（整个立项决策里最关键的一条）

一个机制不可知的ADC，在PFS显著改善的情况下，因为OS(HR 0.98，几乎零获益)未达
统计学显著性，被公司主动撤回已经提交的BLA。这说明在这个适应症里：
1. 光有ORR/PFS不够，监管和临床实践最终看的是OS；
2. 病人在后线的基线状况(多线治疗后)让OS终点天然很难打出显著差异——
   这是任何后线方案(包括本仓库三条路线里的任何一条)都要面对的共同风险，
   不是某个具体机制的问题。

## 二、四代EGFR TKI(路线C)竞品核查

| 药物 | 公司 | 机制 | 关键数据(已核实) |
|---|---|---|---|
| silevertinib(BDTX-1535) | Black Diamond Therapeutics | **共价**，绕开C797的新位点(非经典突变+C797S广谱覆盖) | Fast Track认证；ASCO 2026口头报告C797S队列数据(42例，100mg/200mg剂量组)；**公开物化性质：MW 561.05，TPSA 81.56 Ų，XLogP 0.77**；**CNS Kp,uu 大鼠0.48-0.8/犬0.48(两个文献来源数值不完全一致，已标注)**；在颅内PDX模型中有效 |
| tigozertinib(BLU-945) | Blueprint Medicines | **可逆/非共价**，wild-type-sparing | 对T790M/C797S双突变和三突变均有活性，SYMPHONY临床试验(108例重度预处理患者)，WT抑制相关毒性<10% |
| BBT-176 | Boryung | 可逆ATP竞争 | 1/2期临床(NCT04820023)，对单/双/三突变Ba/F3模型均有活性 |

**重要更正**：早期讨论稿把"主攻可逆正构(reversible_orthosteric)"的precedent
错误归到了silevertinib头上——核查后确认**silevertinib实际是共价药**(只是换了
个不依赖C797的新共价位点)，真正验证"可逆机制能打C797S"的precedent是
**tigozertinib(BLU-945)**。这个区分不是文字游戏：如果route_c要主攻可逆正构
路线，应该对标tigozertinib的SAR和临床数据，而不是silevertinib的。

**silevertinib的真正竞争意义**：它公开披露了真实的脑渗透数据，证明"四代EGFR TKI
+CNS穿透"这个本该是技术壁垒的组合已经被实现并且进入了后期临床——这是本仓库
路线C(`routes/route_c_4th_gen_tki`)的TPP校准基准，详见
[doc/routes/route-c-4th-gen-tki.md](routes/route-c-4th-gen-tki.md)。

## 三、SHP2/SOS1变构抑制剂(路线A)核查

### SHP2(PTPN11)——主攻靶点

- **真实结构**：PDB **5EHR**(SHP099/SHP2复合物)——SHP099结合在由N-SH2/C-SH2/PTP
  三域组成的"tunnel"变构位点，稳定SHP2的自抑制构象；关键相互作用残基
  **Phe113、Arg111、Glu250**。另有PDB **7XBQ**(共价片段图谱)可供参考。
- **临床precedent**：TNO155(Novartis)、RMC-4630(Revolution Medicines)均作用于
  同一tunnel位点。**但领域仍处早期**——TNO155在118例重度预处理患者的单药剂量
  爬坡中，仅20%达到疾病稳定(中位维持4.9个月)，联用策略的真实疗效尚未被充分
  证实。这比"已验证靶点、方法学成熟"这种表述要更谨慎。

### SOS1——备选靶点

- **真实结构**：PDB **6SCM**(BI-3406复合物)——结合在紧邻RAS^cat对接界面的浅
  疏水口袋，喹唑啉环和His905发生π堆积；另有PDB **1NVU**(催化位点参考)、
  **7AVI**(化合物2复合物)。
- **一个有趣的历史细节**：BI-3406的前体化合物BI-68BS最初是EGFR项目的产物，
  后来才被重新定位到SOS1——这类"跨靶点化学型复用"在药化史上并不罕见。

## 四、EGFR降解剂(路线B)核查

- **真实precedent**：CFT8919(C4 Therapeutics，EGFR L858R选择性BiDAC降解剂)，
  对L858R-C797S/L858R-T790M/L858R-T790M-C797S均保留活性。**但比预期更早期**：
  IND刚获批，中国(Betta Pharmaceuticals合作)的1期试验截至核查时"尚未开始入组"。
  这意味着整个EGFR降解剂赛道几乎没有临床验证数据——机会和风险并存。
- **E3连接酶选择**：文献里CRBN(来那度胺/泊马度胺衍生配体)和VHL(如VH032类配体)
  都是常见选择。真实公开案例：MS39(VHL招募，gefitinib骨架)、MS154(首个CRBN
  招募的EGFR降解剂)、SIAIS125(CRBN招募，canertinib骨架+泊马度胺)。
- **CNS渗透是真实的、文献公认的短板**：有综述明确指出"大多数EGFR降解剂体外
  活性良好，但PK性质(尤其是CNS穿透)往往跟不上"——这不是本仓库编的风险描述，
  是这个领域公认的技术瓶颈，也是`routes/route_b_degrader`把CNS Kp,uu提升为
  一级TPP指标的直接依据。
- **来那度胺/泊马度胺对CRBN的真实结合常数**：来那度胺Ki≈3.1 μM，泊马度胺
  Ki≈0.8 μM(不同assay格式下数值有差异，属于该领域常态)。已用于
  `core/ternary_complex.py::evaluate_hook_effect_risk()`的真实数据演示，
  见`doc/routes/route-b-degrader.md`。

## 四.5、EGFR C797S真实晶体结构格局核查(比预想的更稀缺)

- **del19/C797S(路线C primary genotype)没有任何晶体结构**——文献明确提到
  del19/T790M/C797S三突变体"蛋白都没能纯化到结晶级别"。
- **L858R/C797S双突变体(无T790M，路线C另一个primary genotype)也检索不到
  任何PDB条目**。
- 现有全部C797S相关晶体结构(**6LUD**/6LUB/9D3V/9D3W/9XU9/8WD4等)都是
  L858R/T790M/C797S三突变体——也就是路线C的secondary genotype。
  本仓库已下载**6LUD**(L858R/T790M/C797S + 奥希替尼本身的真实共晶，配体YY3
  的37个重原子和奥希替尼C28H33N7O2精确吻合)，存在
  `routes/route_c_4th_gen_tki/structures/`，并登记进数据库`structure_ensemble`表。
- **路线A的SHP2/SOS1真实结构可以直接下载**：PDB 5EHR(SHP099/SHP2 tunnel
  位点)、6SCM(BI-3406/SOS1)都已下载+清洗，配体HETATM(5OD/L7H)确认存在。
  这是路线A相对路线C的一个具体优势——不需要先走计算突变这一步。

## 四.6、拿真实上市/临床药物跑本仓库自己的ADMET流水线(验证计算工具本身可信度)

用RDKit核对过SMILES与PubChem公开分子式精确匹配(奥希替尼C28H33N7O2/499.6，
**tigozertinib(BLU-945) C28H37FN6O3S/556.7**，**silevertinib(BDTX-1535)
C30H30ClFN6O2/561.05**)，跑了一遍`core/admet.py`/`core/synthesis.py`：

| 化合物 | MW | CNS_MPO | SA score | hERG代理 |
|---|---|---|---|---|
| osimertinib | 499.6 | 2.25 | 2.92 | high |
| tigozertinib(BLU-945) | 556.7 | 2.99 | 4.53 | **low** |
| silevertinib(BDTX-1535) | 561.1 | 2.35 | 4.14 | high |

意外但可解释的发现：tigozertinib是唯一hERG代理显示"low"的，因为它的哌啶
氮直接连在嘧啶环上(芳香胺，不是经典烷基叔胺)，碱性被电子不足的杂环削弱，
`core/admet.py`的`has_basic_aliphatic_amine`判断正确识别了这个差异——
这和本仓库`pocket_regions.yaml`"solvent_exposed区域富集弱碱性替代物"的
设计原则吻合，说明这条启发式规则抓住了真实的药物化学模式，不是凭空编的。
详细解读(包括"不要被自己demo分子的CNS_MPO数字误导"这条重要警示)见
`doc/routes/route-c-4th-gen-tki.md`。

## 五、三路线立项逻辑小结

碎片化的战场上，赢的策略是"收敛"而不是"特异"——不要去抢某一个7%的分片，
找能同时覆盖多个分片的收敛点：

- **路线A(SHP2/SOS1)**：收敛点是下游信号节点——MET扩增/HER2扩增/PIK3CA/
  RAS通路重新激活，合计覆盖约25-30%的耐药病人，且和EGFR TKI联用有明确生物学
  逻辑。变构位点结构清楚、共晶多，是三条路线里FEP最友好的一个。
- **路线B(EGFR降解剂)**：收敛点是作用机制——不挑具体突变类型，理论上能同时
  覆盖C797S和其他靶上突变，还能绕开"压过细胞内毫摩尔级ATP"这个可逆/共价抑制剂
  共同的难题。代价是CNS渗透这个真实短板，以及三元复合物预测这个尚不成熟的
  计算化学前沿问题。
- **路线C(四代EGFR TKI)**：技术风险最低、工程实现最完整，但核查后发现是三条
  路线里竞争最激烈的一个——三个真实项目同时在跑，其中silevertinib已经公开
  披露了本该是护城河的CNS渗透数据。保留这条路线的价值在于：它是目前唯一
  有完整端到端计算管线实现的路线(`routes/route_c_4th_gen_tki`)，可以作为
  验证整套方法论(标准化/枚举/打分漏斗/ADMET/决策引擎/DMTA回传)的参照实现，
  路线A/B复用的正是这套共享引擎(`core/`)。

一个贯穿三条路线的原则：**从第一天起就按"联用"设计，不要按"单药"设计**。
一线奥希替尼耐药后已经不存在单药能赢的局面(amivantamab+化疗是现在的标准)。
三条路线的`target_profile.yaml`都新增了CYP药物相互作用风险、与常用联用药物
不重叠的毒性谱、给药方案兼容性这三项TPP指标，就是这个原则的具体落地。
