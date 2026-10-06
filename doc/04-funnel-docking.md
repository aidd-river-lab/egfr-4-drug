# 环节4：多级打分漏斗 L0-L4

代码：[`core/docking.py`](../core/docking.py) [`core/md_stability.py`](../core/md_stability.py) [`core/mmgbsa.py`](../core/mmgbsa.py) [`core/fep.py`](../core/fep.py) [`core/covalent.py`](../core/covalent.py) · 配置：[`config/pipeline.yaml`](../config/pipeline.yaml) · **状态：L1对接2026-10起真实可跑(见`scripts/setup_docking_env.sh`)；L2质控逻辑真实；L3短程小规模MD 2026-10起真实可跑(配体结合复合物，10ps demo)，真实项目量级(15ns×3副本)/L4(FEP)仍是诚实占位，需要GPU**

## 核心思想：每升一级精度提高一个量级，通量降低两个量级

原设计文档开篇就警告"用Vina做最终排序，等于用卷尺量头发丝"——快速对接的打分
函数误差在2-2.5 kcal/mol量级，而真正区分两个相似分子活性差异所需的精度是
零点几kcal/mol。所以漏斗设计的关键不是"每一级都要很准"，而是"前面几级只做
排除，只有足够贵的后面几级才用来排名"。

| 层级 | 方法 | 工具 | 通量/天 | 误差 | 作用 | 代码 |
|---|---|---|---|---|---|---|
| L0 | 药效团+形状+属性过滤 | ROCS/Pharmit/RDKit | ~10^8 | — | 从巨型库粗排除 | `core/enumerate.py`+`core/admet.py` |
| L1 | 快速对接 | AutoDock Vina/Gnina | ~5×10^5 | 2.5 | 砍掉明显不贴合的，**不排名** | `core/docking.py` |
| L2 | 精细对接+重打分 | Glide XP/Gnina CNN/共价对接 | ~10^4 | 2.0 | 生成可信结合构象供L3用 | `core/docking.py`/`core/covalent.py` |
| L3 | MM-GBSA+短程MD稳定性 | OpenMM/GROMACS | ~300 | 1.3 | 排除"对接分高但不稳定"的假阳性，**性价比最高的一步** | `core/md_stability.py`/`core/mmgbsa.py` |
| L4 | FEP相对结合自由能 | OpenFE/FEP+/Amber TI | ~30/GPU节点 | 0.75 | 同系列SAR的真正预测，决定合成哪个 | `core/fep.py` |

完整配置（含每层的pass_criteria、工具列表）见 `config/pipeline.yaml`。

## L1-L2：对接（诚实占位，含一个修复过的真实bug）

`meeko`（配体/受体PDBQT准备）需要 Python≥3.10（本环境3.9.6，meeko内部用了
`match`语句），`vina` 需要系统预装Boost C++库——两者本环境都装不上。
`core/docking.py` 用lazy import包起来，import模块本身不会报错，真正调用时才
诚实返回 `ok=False` + 具体缺什么。

**这个模块存在的第二个目的是修正姊妹项目 `egfr-pipline` 的一个真实bug**：
那个项目的"受体PDBQT准备"函数实际上就是 `cp receptor.pdb receptor.pdbqt`——
只换了文件后缀，完全没做AutoDock需要的原子类型/部分电荷/可旋转键计算，
导致那个项目产出的全部对接 `binding_affinity` 都是 `0.0`（不是真实Vina打分，
是PDBQT格式不合法时Vina的某种退化行为）。`core/docking.py::prepare_receptor_pdbqt()`
改用meeko官方推荐的 `mk_prepare_receptor.py` CLI 做真正的转换。

```python
from core.docking import prepare_receptor_pdbqt, prepare_ligand_pdbqt, run_vina_docking

prepare_receptor_pdbqt("receptor_clean.pdb", "receptor.pdbqt")
# ok=False: "mk_prepare_receptor.py 不在PATH上(meeko未安装...)"

run_vina_docking(receptor_pdbqt=..., ligand_pdbqt=..., center=(-14.2, 33.5, 22.8), box_size=(20,20,20), out_pdbqt=...)
# ok=False: "vina 未安装或不兼容当前Python版本...（还需要系统预装Boost库）"
```

要在新环境接入：`bash scripts/setup_docking_env.sh` 一键搭建`.venv310`(Python 3.10 +
meeko + vina)。这个脚本不是简单的pip install——记录了三个真实踩过的坑(vina的
setup.py硬编码的Boost查找路径、老版本C++标准和新版Boost的类型别名冲突、缺swig)，
详见脚本内注释和`core/docking.py`模块docstring。

## 2026-10：真实对接已经跑通，不再是纯占位

用这套环境对真实下载的PDB结构(路线A的5EHR，路线C的6LUD)做了端到端真实对接：
受体用`mk_prepare_receptor.py`正规转换(处理了真实晶体结构的altloc问题)，配体
用RDKit生成3D构象+MMFF优化再转PDBQT，口袋中心取自真实共晶配体的坐标质心
(不是猜的)。路线A全部16个候选、路线C的12个决策批次候选都跑出了真实、有区分度
的结合能分数(范围大约-7.5到-10.7 kcal/mol，没有一个是0.0)。

**一个有效性检验**：把奥希替尼自己的SMILES重新对接进6LUD(它自己的共晶结构)，
得到-7.76 kcal/mol，和其它候选分子没有显著差异——这符合已知药理学：奥希替尼
真正的高效力来自和Cys797形成的共价键，而C797S恰好去掉了这个共价靶点，Vina只能
打出它的非共价结合姿势分数，自然不会特别突出。这是一个合理性检验，不是精确的
构象重现验证(没有做对接姿势和晶体姿势的RMSD比对)。

**一个比预期更重要的发现**：把真实对接分数和环节8算出来的desirability(纯ADMET，
不含任何结合信息)做Spearman相关，路线C的12个候选相关系数只有0.25(p=0.43，
不显著)——**ADMET打分和真实结合强度几乎没有关系**，这正是本仓库从一开始就
强调"不能只靠环节6/8做最终决策，必须有环节4的打分漏斗"的原因，现在有真实数据
支撑这句话了。两条路线里，`select_batch()`的control档(从"desirability较低"的
池子里随机抽的)都抽到了真实对接分数名列前茅的分子(路线A的RTA-0005综合分数
全场最低但对接分数-10.74全场最佳；路线C的RTC-0034/RTC-0134同样排进前五)——
这是"假阴性复活"机制在真实数据上的具体验证，不是假设性的设计理念。

## L2.5：共价对接（chemistry_route=covalent_new_site/covalent_pan_mutant_broad专用）

仅当 `target_profile.yaml` 的 `chemistry_route.covalent_warhead_enabled=true`
时启用。共价对接需要Schrödinger CovDock/AutoDock4-covalent/Rosetta，本环境
未安装，`run_covalent_docking_stub()` 诚实占位。

**真实可跑、经过测试的部分**是Bürgi-Dunitz攻击角判据——弹头碳原子攻击亲核
残基时，几何上要求攻击角接近105°（纯几何判断，不需要跑过共价对接就能验证）：

```python
from core.covalent import evaluate_attack_geometry, CovalentDockingResult

result = CovalentDockingResult(
    ligand_id="DEMO", nucleophile_residue="Lys745",
    warhead_carbon_to_nucleophile_distance_angstrom=3.2, attack_angle_degree=103,
    ki_nm=50, kinact_per_second=0.02,
)
evaluate_attack_geometry(result)
# {"geometry_ok": True, "kinact_over_ki": 0.0004, "note": "共价药的IC50没有意义；排序/决策请用kinact_over_ki..."}
```

**关键认知**（原设计文档"环节4.3末尾"）：共价药的IC50没有意义，要看
`kinact/KI`。`KI`(可逆亲和力)决定选择性，`kinact`(成键速率)决定共价效率，
两者优化方向不同（提高KI靠结合位姿，提高kinact靠弹头几何对齐和亲核体pKa），
所以 `CovalentDockingResult` 把两者分开存，不合并成一个数。

## L3：MM-GBSA + MD稳定性

真实项目需要的量级(3副本×20ns)需要GPU，本环境没有，`core/mmgbsa.py::
run_mmgbsa_stub()` 诚实占位。

**2026-10更新：短程、小规模的真实MD现在能跑了**——`openmm`+`pdbfixer`纯pip
可装(`.venv310`)；配体结合复合物还需要`openff-toolkit`+`openmmforcefields`+
真实AmberTools(antechamber二进制)，这几个pip装不了(openff-toolkit需要
Python≥3.11)，装在一个独立的conda环境里(`mamba install -c conda-forge
ambertools openff-toolkit openmm openmmforcefields pdbfixer`，详细步骤见
`scripts/setup_docking_env.sh`第7-8步)。

真实跑过一次：`core/md_stability.py::run_protein_ligand_complex_md()`，
奥希替尼(**真实6LUD晶体坐标**，不是对接预测的姿态——用
`AllChem.AssignBondOrdersFromTemplate()`把已验证过的真实SMILES的键级信息
转移到真实晶体坐标上，比用对接姿态更贴近真实结合模式) + 6LUD受体(C797S
三重突变)，5051原子，真实能量最小化(526774→-40599 kJ/mol)，真实10ps轨迹，
配体RMSD轨迹`[1.25, 1.27, 1.36, 1.46, 1.45, 1.92, 1.6, 1.75, 1.86, 1.5]`，
均值1.54Å，没有发散——**这个真实晶体姿态在短程MD下是稳定的**，和"奥希替尼
确实能非共价结合C797S突变体，只是结合力不如共价焊接牢"这个已知生物学事实
一致，没有出现假阳性迹象。

**这证明了工具链本身是通的，不代表"L3已经真实可用"**：这只是10ps的demo，
真实项目要看的是15ns窗口的均值(pass criteria下面写的2.5Å阈值)，量级差了
三个数量级；`hinge_hbond_occupancy_pct`/`target_anchor_occupancy_pct`这两个
字段的真实计算逻辑还没实现(需要按残基名追踪氢键距离，目前诚实留空，不编造
数字)；电荷方案用的是gasteiger(RDKit内置，免量子化学)，不是生产级AM1-BCC。
按`run_protein_equilibration_md()`实测的CPU吞吐量推算，15ns×3副本这个量级
需要数十小时到几天，这是`run_md_stability_stub()`继续保留"需要GPU"占位的
真实依据，不是没去试。

Pass criteria（`config/pipeline.yaml` `funnel.L3.pass_criteria`）：
- 配体RMSD(最后15ns均值) < 2.5 Å
- 铰链氢键占有率 > 60%
- target_anchor(C797S路线里指Ser797-OG)占有率 > 40%——**注意这不是硬门槛**

### 一个"致命坑"：Ser797氢键不是硬门槛

`config/pipeline.yaml` 的 `ser797_hbond_policy`：

```yaml
ser797_hbond_policy:
  is_hard_gate: false
  max_score_weight_pct: 15
  must_use_md_occupancy: true   # 不能只看单帧对接位姿
```

C797S四代药的活性主要来自对整个ATP口袋的高亲和力贴合，与Ser797直接氢键是
加分项不是必要条件。把它设成硬门槛会误杀好分子——这是原设计文档"致命坑清单"
的第10条，之所以单独强调，是因为直觉上"这个突变体缺了一个关键氢键供体"很容易
被想当然地做成一票否决的门槛，但实际药化数据不支持这个假设。

## L4：FEP相对结合自由能

需要OpenFE/OpenMM或商业FEP+，单对微扰2-8 GPU小时，50对的map需要
200-400 GPU小时，本环境无GPU。`core/fep.py::run_fep_stub()` 诚实占位。

**真实可跑、经过测试的部分**是质控逻辑——`FEPMapResult.qc_pass()` 实现了
原设计文档"环节4.4"的硬规则：

```python
from core.fep import FEPMapResult, FEPPerturbation

bad_map = FEPMapResult(genotype="del19_C797S", perturbations=[
    FEPPerturbation("A", "B", predicted_ddg_kcal_mol=2.5, uncertainty_kcal_mol=0.3,
                    experimental_ddg_kcal_mol=-1.0),
])
bad_map.qc_pass()
# {"ok": False, "reason": "MUE=3.50 kcal/mol 超过阈值1.5，整张map作废，回去查位姿或质子化态"}
```

规则：(1) 没有实验锚点的map不能用来做决策；(2) 预测vs实测的MUE超过
`config/pipeline.yaml` 的 `funnel.L4.qc.mue_reject_threshold_kcal_mol`(1.5
kcal/mol)，整张map作废——不是挑几个好的点用，是整张map都不可信，因为MUE超标
通常意味着受体结构/质子化态系统性错了，局部修补没有意义。
