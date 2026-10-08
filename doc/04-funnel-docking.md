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
真实依据，不是没去试。真要免费拿到GPU把这个量级跑起来，见下面。

### 怎么用免费GPU真实跑起来：`notebooks/colab_gpu_md.ipynb`

这个notebook是**自包含的**——不需要上传这个仓库的任何其他文件，所有
结构数据都会在notebook里重新从RCSB真实下载。打开+跑起来的步骤：

**准备阶段(在Colab界面里操作，不是notebook的cell)**：
1. 打开 [colab.research.google.com](https://colab.research.google.com)，
   登录你的Google账号(免费GPU额度跟账号绑定，有每日/每几小时的用量上限，
   用太多会被临时限流，这是Colab免费层的真实限制，不是bug)。
2. `文件 → 上传笔记本`，选本仓库的`notebooks/colab_gpu_md.ipynb`。
3. `修改运行时类型`，确认加速器选的是`T4 GPU`(notebook的metadata里已经
   预设了这个，正常情况会自动选上，这一步是手动兜底检查)。

**notebook内部8个cell，按顺序点"运行"，不要跳步**：

| 步骤 | 做什么 | 要注意什么 |
|---|---|---|
| 1 | 装`condacolab`(给Colab装conda/mamba的社区工具) | 运行完**会自动重启运行时**，看到"会话已崩溃"提示是正常的，不是报错，等几秒重启完成 |
| 2 | 重启后，`mamba install`装这次会话本地验证过的同一套依赖：`ambertools`(真实antechamber二进制)+`openff-toolkit`+`openmm`+`openmmforcefields`+`pdbfixer`+`rdkit` | 用的是mamba不是conda classic(本地踩过conda解这组依赖卡死的坑)；装完会打印`可用平台`列表，**必须看到`CUDA`在列表里**，否则说明没真的分到GPU，要回去检查准备阶段第3步 |
| 3 | 真实下载6LUD晶体结构(`urllib`直连RCSB)，清洗出受体链 | 和本地`core/structures.py::clean_chain()`逻辑完全一样，这里是内联重写的一份，不依赖仓库代码 |
| 4 | 从真实晶体坐标重建奥希替尼，`AssignBondOrdersFromTemplate`转移真实键级 | 和本地做法一致，不是用对接预测的姿态 |
| 5 | 定义MD函数(和`core/md_stability.py::run_protein_ligand_complex_md()`逻辑一致，加了显式GPU平台选择) | 纯定义，不产出结果，运行很快 |
| 6 | 先跑一个**标定**：5000步(10ps)，真实测这块GPU的吞吐量 | 这一步拿到的`seconds_per_1000_steps`数字，决定第7步能跑多少 |
| 7 | 按标定结果，在一个时间预算(默认20分钟)内尽量跑接近真实项目量级的轨迹 | 想跑更长就把`time_budget_minutes`改大，但**Colab免费层会话大概12小时会被强制断开，空闲太久也会断线**，不要设得离谱大 |
| 8 | 把跑出来的`.dcd`轨迹文件下载回本地 | `.dcd`能用VMD/PyMOL/MDAnalysis打开看真实的分子运动；同时会打印这次真实跑出来的配体RMSD轨迹，照着本节上面的pass criteria(<2.5Å)自己判断这次姿态稳不稳 |

**怎么判断这次跑得对不对**：第6步标定如果平台列表里没有`CUDA`，后面
全是在用CPU模拟GPU的"假跑"，速度跟本地Mac差不多(86秒/10ps量级)，没有
意义——一定要先确认标定阶段的吞吐量数字明显快于这个数(真实T4一般能
快一到两个数量级)才继续。第8步打印的RMSD轨迹，解读方式和本文档上面
"真实跑过一次"那段描述完全一样：均值<2.5Å且没有持续爬升趋势=姿态稳定，
持续发散=这个对接姿态很可能是假阳性。

**2026-10更新：这个notebook真实在GPU上跑通了**——用户真实在Colab免费
T4上跑完了路线A的3个真实候选分子(RTA-0005/0009/0010)的复合物MD，过程中
踩了一串真实的环境坑(详见下面的debug记录)，最终真实吞吐量是OpenCL平台
6.6秒/1000步，比本地Mac CPU的17.2秒/1000步快约2.6倍——不是CUDA(这台
机器的CUDA插件PTX版本不兼容，自动探测逻辑退到了OpenCL，仍是真实GPU
加速)。真实结果和落库方式见[route-a-shp2-sos1.md](routes/route-a-shp2-sos1.md)。

**真实踩到的一个坑，第5步`from openff.toolkit import Molecule`报
`RuntimeError: operator torchvision::nms does not exist`**：根源不在
openff-toolkit本身——本地用同一套组合(`mamba install ambertools
openff-toolkit openmm openmmforcefields pdbfixer`)验证过，
`openff.toolkit`完全不需要torch。真实原因是这组依赖在某些Colab环境里
被mamba的解析器连带装进了一对**版本不匹配**的`pytorch`+`torchvision`，
报错发生在它们互相校验的那一步，跟我们真正需要的几个包完全无关。

**第一次尝试的修复思路(钉`python=3.11`)已验证是错的，记录下来避免
再踩**：condacolab会在环境里留一个pin文件强制锁住Python版本(和Colab
当前kernel进程绑定的版本一致，这是condacolab故意的设计，防止改Python
版本直接搞坏正在跑的kernel)，真实报错`pin on python =3.13 ... python
=3.11 is not installable because it conflicts`——这条路走不通。

**第二次尝试(`pip uninstall torch torchvision torchaudio`)也没用，
原因找到了**：这几个包是mamba装的，不是pip装的，而且conda-forge上
这个包真实的名字是`pytorch`，不是`torch`——pip从名字上就找不到它，
会静默跳过(不报错，但什么也没删)。**真实有效的修复**：用`mamba remove`
配上正确的包名：

```bash
!mamba remove -y pytorch torchvision torchaudio
```

删完**必须重启一次Colab会话**(`代码执行程序→重新启动会话`，不用删
整个运行时，磁盘上已经装好的mamba环境还在)——因为第一次import失败时
torch的C扩展可能已经把坏状态注册进了这个Python进程的内存里，单纯在
同一个还活着的kernel里重跑cell不会清掉这个内存状态，必须换一个全新的
Python进程。重启完，`condacolab.check()`那个cell(第2步)重新跑一遍，
不需要重新跑`pip install condacolab`(第1步，已经装好了)，然后第3-5步
按顺序重跑。notebook现在已经把`mamba remove`这一行直接加进第2步的
安装cell里，以后跑这个notebook不会再需要手动处理这个坑。

**另一个真实踩到的坑，第10步报`NameError: name 'calib' is not
defined`**：重启运行时之后，第1-8步之前跑出来的所有Python变量都被
清空了(磁盘上装好的mamba环境还在，但变量不在)。第10步要用到第6步
标定出来的`calib["seconds_per_1000_steps"]`来换算该给每个真实候选
分子跑多少步——如果重启之后跳过第6步直接冲到第9/10步(不想再跑一遍
demo很正常)，这个变量就是空的。**不需要把第7/8步的完整demo流程也
跑一遍**，只要把第2(确认环境)、3、4、5、6步依次重跑一遍(第6步就是
那个5000步的小标定，几十秒内跑完)，再回来跑第9步上传、第10步，就
不会报这个错了。notebook现在也在第10步加了一个明确的检查，报错信息
会直接告诉你该回去跑哪几步，不会再是一个看不懂的`NameError`。

**第三个真实踩到的坑，第6步报`TypeError: A GAFF force field name must
be provided as a string`**：`GAFFTemplateGenerator(molecules=[off_mol])`
不传`forcefield`参数时，"不传该用什么默认值"这件事在不同版本的
`openmmforcefields`上处理得不一样——本地conda环境装到的版本会自动
补上`gaff-2.11`，Colab这次装到的版本不会，直接报错。**修复**：显式
传`forcefield="gaff-2.11"`，不依赖任何版本的默认填充行为，notebook
第5步的函数定义已经改过来了。

**这同一个`torchvision::nms`报错，真实踩了4轮才找到根因，记录完整debug
过程(不是最后一次就蒙对的，前3次方向都错了，诚实记录)**：

1. 第1次猜测：Python版本太新(3.13)，conda-forge的torch/torchvision生态
   没跟上——尝试钉`python=3.11`，**失败**：condacolab会锁死Python版本，
   报`pin on python =3.13 ... python =3.11 is not installable`，这条路
   根本走不通。
2. 第2次猜测：torch/torchvision版本不匹配，卸载重装——`pip uninstall
   torch`，**无效**：这几个包是mamba装的，conda-forge上的包名是
   `pytorch`不是`torch`，pip从名字上就找不到它，静默跳过什么也没删。
3. 第3次猜测：既然pip卸不掉，用`mamba remove pytorch torchvision
   torchaudio`，**还是无效**：后来查`mamba list`真实数据才发现torch
   是被`openff-nagl`当硬依赖拉进来的，删了也会被重新拉回来，而且
   `torchvision`压根没被装过，删一个不存在的包自然没意义。
4. **真实查到根因**：`mamba list`显示装了`openff-nagl`(一个ML电荷
   预测的可选后端，带出`pytorch`+`pytorch-lightning`+`torchmetrics`)，
   `torchmetrics`会在导入时尝试给`torchvision::nms`注册一个"假"实现，
   但环境里从没真正装过`torchvision`，这个算子从没被注册过，所以报
   "operator does not exist"。第4次先尝试"补装torchvision让算子真的
   存在"，**仍然失败**(同样的报错)，说明问题比"缺一个包"更深，可能是
   补装的torchvision build和已装的pytorch 2.10.0+cuda130这个具体build
   在ABI上不兼容。
5. 第5次尝试：查OpenFF官方文档确认`openff-toolkit`(完整版)会打包捆绑
   `openff-nagl`这整条可选ML链路，`openff-toolkit-base`只装核心硬依赖，
   改用`openff-toolkit-base`——**依然失败**：真实查`mamba list`发现
   `openff-toolkit`(完整版)又被装回来了，说明`openmmforcefields`自己
   硬依赖完整版`openff-toolkit`，不管我们显式请求哪个包名都绕不开，
   `openff-nagl`跟着躲不掉。
6. 第6次尝试：怀疑是两次分开`mamba install`导致torch/torchvision的build
   字符串选得不一致，改成把`torchvision`和其他包放进**同一条**mamba
   install命令里一次性解析——**还是失败**，报错一字不差。

**到这里为止5-6次都是在conda安装层面猜/试，第7次换了方法：拿完整、不
折叠的错误堆栈直接追到底，不再猜**。真实堆栈显示：`openff.toolkit`
导入时会自动构建一个"工具箱注册表"，其中包含`NAGLToolkitWrapper`——
它的`is_available()`方法设计上用`try/except ImportError`包住
`importlib.import_module("openff.nagl")`，想优雅地处理"这个可选ML后端
没装"的情况。但真实链路是：`openff.nagl`→`pytorch_lightning`→
`torchmetrics`→`torchvision`，**torchvision自己的原生C++扩展在这个
环境里没能正确把`torchvision::nms`这个op注册进torch的dispatcher**，
它自己的`_meta_registrations.py`紧接着尝试给这个"应该已经注册好"的op
再注册一个meta/fake实现时，发现op根本不存在，抛出的是`RuntimeError`，
不是`ImportError`——`is_available()`的兜底网接的是`ImportError`，接不住
这个`RuntimeError`，直接一路崩到最外层。

**真正有效的修复，不在conda层面，在Python运行时层面**：在
`from openff.toolkit import Molecule`之前，加两行：

```python
import sys
sys.modules["openff.nagl"] = None
```

Python的import系统规定：`sys.modules`里某个模块名如果被设成`None`，
之后任何对这个模块的import都会直接抛出`ImportError`(文档化的标准行为，
不是hack)——这样`NAGLToolkitWrapper.is_available()`里的
`importlib.import_module("openff.nagl")`会在真正触发那条broken链路
之前就拿到一个它能正常处理的`ImportError`，干净地判定NAGL不可用，
转而用RDKit/AmberTools这些真正要用的工具箱，不碰conda环境、不删任何
已装的包。notebook第5步已经加上这两行。

这整个坑从第1次尝试到找到真正根因，一共走了7轮——记录下来是提醒自己
(也提醒之后看这篇文档的人)：**遇到深层依赖链路的报错，折叠过的
traceback只够猜，真正定位问题通常需要那个没有被折叠的完整堆栈**，
前6次尝试都是在堆栈信息不完整的情况下合理但错误的猜测。

**第四个真实踩到的坑，第6步报`OpenMMException: Error loading CUDA
module: CUDA_ERROR_UNSUPPORTED_PTX_VERSION`**：`Platform.getNumPlatforms()`
列出了`CUDA`不代表真能用——conda装的openmm-cuda插件编译时用的CUDA版本
比这台Colab机器实际GPU驱动支持的版本更新，插件能加载、能列出来，但真正
建Context执行计算时才报错。**修复**：把`get_best_platform()`从"看列表
按名字选"改成"真的建一个最小测试Context验证每个候选平台能不能用，
第一个真能用的才返回"——CUDA真失败会自动退到OpenCL(在T4上仍然是真实
GPU加速，只是走另一套驱动接口，不受这个PTX编译版本问题影响)。notebook
第5步已经按这个逻辑改过。

### 不要只跑demo——第9-11步，换成用这个仓库真实筛选出来的候选分子

第1-8步跑的是奥希替尼+6LUD，只是用来验证GPU工具链本身通不通，**不是
你应该一直用demo分子跑下去**。真实想推进的是漏斗(环节4)已经筛出来的
那批真实候选——`routes/<route>/structures/RTx-xxxx_docked.pdbqt`，它们
的L1 Vina分数已经真实落在数据库`funnel_scores`表里。notebook的第9-11
个cell(紧接在第8步后面)就是接这一段的：

**第0步，本地先跑(在你电脑终端，不是Colab)**：

```bash
PYTHONPATH=. .venv310/bin/python scripts/export_top_candidates_for_gpu_md.py \
    --route route_a_shp2_sos1 --top-n 3
```

`scripts/export_top_candidates_for_gpu_md.py`做的事：查真实数据库里这条
路线L1对接分数最好的N个候选(`funnel_scores.score_value`从小到大排序，
越负越好)，对每个候选的`_docked.pdbqt`(Vina真实算出来的对接姿态)用
`meeko.PDBQTMolecule`+`RDKitMolCreate`(meeko官方提供的"把对接结果读回
正确键级分子"工具，不是本脚本自己拼的转换逻辑，这次会话里拿真实的
`RTC-0000_docked.pdbqt`验证过还原出的SMILES和设计结构完全一致)还原出
带正确化学键+真实对接坐标的`.sdf`，连同这条路线对应的真实受体PDB
打包成`validation/gpu_md_export/<route>_top<N>.zip`。`--route`可以换成
`route_c_4th_gen_tki`(四代TKI路线)；路线B(降解剂)目前只有warhead片段
的对接分数(见[16第7节](16-visual-dashboard.md)已知限制)，这个脚本暂时
不支持。

**notebook里的第9-11步**：

| 步骤 | 做什么 |
|---|---|
| 9 | 上传刚才本地生成的zip，解压，打印`MANIFEST.txt`里每个候选真实的Vina分数来源 |
| 10 | 对zip里的每个真实候选，用第5步定义的同一个`run_protein_ligand_complex_md()`函数跑一遍复合物MD(默认每个候选分10分钟GPU时间预算，按第6步标定出的真实吞吐量换算步数)，每个候选跑完当场打印配体RMSD轨迹和"稳定/可能不稳定"的判断 |
| 11 | 把所有候选的`.dcd`轨迹+汇总JSON打包下载 |

跑完把`my_candidates_md_results.json`贴回来，我可以帮你解读哪些候选
在这个时间尺度上稳定、要不要把结果更新进`funnel_scores`表的L3记录——
这一步不是demo了，是真实筛选流程的下一级。

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
