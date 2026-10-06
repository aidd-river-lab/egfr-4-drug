# 环节7：合成可及性

代码：[`core/synthesis.py`](../core/synthesis.py) · **状态：第一层(SA score)真实可跑；第二层(逆合成)2026-10起真实可跑(AiZynthFinder，`.venv310`)；FTO诚实占位(需要真实专利数据库)**

## 三层筛查设计

原设计文档把合成可及性评估分三层，精度和成本依次递增：

1. **启发式分数**（秒级）——RDKit Contrib的SAscore，本模块真实实现。
2. **计算机逆合成**（分钟级，需要训练好的单步反应模型+可购建块库）——AiZynthFinder，
   本模块诚实占位。
3. **药化人工路线规划**（人工评审，不在代码范围内）。

## 第一层：SA score（真实可跑）

```python
from core.synthesis import compute_sa_score

compute_sa_score("COc1cc(N(C)CCN(C)C)c(NC(=O)C=C)cc1Nc1nccc(-c2cn(C)c3ccccc23)n1")  # 奥希替尼
# {"ok": True, "sa_score": 2.92, "interpretation": "1=容易, 10=很难"}
compute_sa_score("CCO")  # 乙醇
# {"ok": True, "sa_score": 1.98, ...}
```

SAscore来自RDKit安装目录下的 `Contrib/SA_Score/sascorer.py`——这不是一个标准
可import的子包，`core/synthesis.py` 顶部手动把这个路径加进 `sys.path` 才能
`import sascorer`（导入顺序不能颠倒，见模块代码注释）。

1(容易合成)到10(很难合成)的量级，只能用来排除明显离谱的候选，**不能单独当
决策依据**——这也是为什么 `workflows/run_discovery_round.py` 把它和CNS MPO
一起喂进 `desirability_score()` 做几何平均，而不是单独用它筛选。

## 第二层：计算机逆合成

AiZynthFinder纯pip可装(`.venv310/bin/pip install aizynthfinder`，不需要像
vina/pyrosetta那样处理Boost/PATH的坑)。真实使用还需要下载官方预训练的
单步反应模型+可购建块库(USPTO历史反应数据训练的模型+约1700万个真实ZINC
可购买分子，合计约1.2GB，不随代码库提交)：

```bash
.venv310/bin/download_public_data validation/aizynth_data
# 会自动生成 validation/aizynth_data/config.yml 指向下载好的全部模型文件
```

**2026-10更新：真实跑过一次，不再是占位**。用这套真实模型+建块库对奥希替尼
(已验证过的真实SMILES)做了一次真实的逆合成搜索：

```python
from core.synthesis import run_retrosynthesis

run_retrosynthesis(osimertinib_smiles, "validation/aizynth_data/config.yml")
# {"ok": True, "is_solved": True, "number_of_steps": 4,
#  "number_of_routes_explored": 100, "number_of_solved_routes": 23, "top_score": 0.975,
#  "precursors_in_stock": "C=CC(=O)Cl, CNCCN(C)C, COc1cc(F)c([N+](=O)[O-])cc1N,
#                          Cn1cc(-c2ccnc(Cl)n2)c2ccccc21",
#  "precursors_not_in_stock": "", "search_time_seconds": 17.5, ...}
```

20秒内搜索了458个节点，找到99条候选路线、24条完全解析到库存建块的路线，
最佳路线4步合成、神经网络打分0.975(满分1.0)。**值得注意的一个细节**：
排名最高的路线里，第一个precursor是`C=CC(=O)Cl`(丙烯酰氯)——这**正是真实
化学里用来安装丙烯酰胺弹头的标准试剂**，AiZynthFinder从历史反应数据里
学出来的这条路线和真实合成化学的直觉吻合，不是瞎编的巧合。

**这是真实结果，但要知道它的真实边界**：这是神经网络从USPTO历史反应数据里
学出来的"看起来合理"的路线，不代表这条路线真的有人验证过能跑通、不代表
产率好，也不代表这是唯一/最优路线——真实项目里这种结果是给合成化学家的
起点建议，不是可以直接照做的操作手册，见`run_retrosynthesis()`函数docstring。

主.venv(Python 3.9)没装aizynthfinder，诚实占位：

```python
run_retrosynthesis_stub(smiles)
# {"ok": False, "reason": "AiZynthFinder 未安装/未配置模型...", "route_found": None}
```

## FTO (Freedom to Operate) 检查

读取 `config/patent_landscape.yaml`——这个文件目前只有占位示例数据
(`status: placeholder_not_real`)，所以 `check_fto()` 现在的返回值永远是
`"unknown_no_real_patent_data"`，**不能当作真实FTO结论**：

```python
check_fto(smiles)
# {"ok": True, "fto_status": "unknown_no_real_patent_data", "note": "patent_landscape.yaml 里还没有真实专利数据..."}
```

真实数据接入后，这里才应该做Markush SMARTS匹配（当前函数体里有这部分逻辑的
占位分支，等真实专利数据接入后生效）。

## 数据模型

`SynthesisRecord` 字段照抄原设计文档"环节7.3"的要求：

```python
@dataclass
class SynthesisRecord:
    compound_id: str
    smiles: str
    sa_score: float | None = None
    n_steps: int | None = None                      # 需要逆合成工具或人工路线规划
    longest_linear_sequence: int | None = None
    key_bb_cas: str | None = None
    bb_lead_time_days: int | None = None
    estimated_cost_per_20mg: float | None = None
    chiral_resolution_needed: bool | None = None
    scale_up_risk: str | None = None                # 自由文本，比如"含叠氮步骤"
    fto_status: str | None = None                    # in_claim / near_claim / clear
```

大多数字段需要逆合成工具或药化人工填写，这个dataclass只定义schema，不是每个
字段都有自动计算方法——没有数据时字段保持 `None`，不为了"看起来完整"去瞎填
默认值。落库表：`db/schema.sql` 的 `synthesis_records`。

## 和环节8的衔接

`target_profile.yaml` 的 `decision_engine.pareto_objectives` 里有一项
`synthesis_n_steps`（对应 `SynthesisRecord.n_steps`，lower_is_better）——
这个指标目前因为第二层逆合成未接入，实际数据永远是 `None`，在真实跑
`core/decide.py` 的Pareto分析前需要先接入AiZynthFinder或人工填写这一列，
否则这个目标维度没有数据可用。
