# 环节6：ADMET描述符 + 结构警示引擎

代码：[`core/admet.py`](../core/admet.py) · 配置：[`config/structural_alerts.yaml`](../config/structural_alerts.yaml) · 测试：[`tests/test_structural_alerts.py`](../tests/test_structural_alerts.py) · **状态：真实可跑，有明确标注的精度边界**

## 精度边界声明（用这个模块做真实决策前必须读）

真正的cLogD7.4需要按pH7.4下的电离态分布算分配系数，这需要一个可靠的pKa预测器
（ChemAxon/Marvin或训练好的pKa模型）。本模块**没有**接入pKa预测，所以：

- `clogp` 是RDKit `Crippen.MolLogP`，这是**中性分子**的辛醇/水分配系数，
  不是cLogD7.4。
- `clogd_approx` 直接复用 `clogp` 的数值作为近似占位，会**系统性高估**带正电荷
  胺类分子在生理pH下的真实分配系数（质子化后极性更强，logD应该更低）。
- `pka_proxy` 不是数值pKa预测，只是"是否含有典型强碱性脂肪叔胺"的子结构标记，
  用来粗略估计hERG/CNS渗透的方向性风险，不能替代真实pKa计算。
- `cns_mpo` 因此也是近似值：MW/TPSA/HBD三项是精确计算，cLogP/cLogD两项用同一个
  MolLogP代入、pKa项用上面的粗代理——这会让分数系统性偏乐观或偏保守，具体方向
  因分子而异。**真实项目决策前必须换成有pKa支持的工具重新算一遍。**

这些边界在模块里每个相关函数的docstring里都重复声明了一遍——宁可啰嗦，
也不要让人在不知情的情况下把这里的数字当成最终结论。

## 描述符计算

```python
from core.admet import compute_descriptors

compute_descriptors("COc1cc(N(C)CCN(C)C)c(NC(=O)C=C)cc1Nc1nccc(-c2cn(C)c3ccccc23)n1")  # 奥希替尼
# {"mw": 499.62, "clogp": 4.51, "tpsa": 87.55, "hbd": 2, "hba": 8, "fsp3": 0.238, ...}
```

MW/TPSA/HBD/HBA/Fsp3都是RDKit精确计算的值，和已知的奥希替尼数据一致
（分子式C28H33N7O2，MW499.62，经过InChIKey交叉验证，见 `doc/01-standardization.md`）。

## CNS MPO（Pfizer六参数复合分）

Wager et al. 2010/2016的desirability function公式，每个参数算0-1分，总分
范围0-6，阈值4.0：

```python
from core.admet import compute_cns_mpo

desc = compute_descriptors(osimertinib_smiles)
compute_cns_mpo(desc)
# {"cns_mpo": 2.25, "pass_threshold_4_0": False, "caveat": "clogd/pka为近似值..."}
```

奥希替尼算出2.25，低于4.0阈值——方向上合理（奥希替尼本身不是为CNS渗透优化的
分子），但由于上面的精度边界，这个具体数值不能当作精确结论。

## 结构警示引擎

每条SMARTS规则都用正/负对照分子真实验证过语法和匹配行为，验证脚本就是
`tests/test_structural_alerts.py`（结构警示库yaml顶部注释里承诺的那份测试，
不是手写猜测的）。9条规则，覆盖的毒性机制和已验证的对照：

| 规则 | 机制 | REJECT/OPTIMIZE | 正对照(应命中) | 负对照(不应命中) |
|---|---|---|---|---|
| `aniline_unmasked` | 代谢活化为醌亚胺的特异质肝毒性 | OPTIMIZE | 苯胺 | 酰胺化苯胺、对/邻位卤代苯胺、2-氨基嘧啶/吡啶类铰链胺 |
| `nitroaromatic` | 还原活化致突变 | REJECT | 硝基苯 | 脂肪族硝基 |
| `furan/thiophene/pyrrole_alpha_unsubstituted` | CYP氧化生成活性中间体 | OPTIMIZE | 未取代呋喃/噻吩/吡咯 | α位取代的同类 |
| `michael_acceptor_offtarget` | 非特异共价/GSH加合 | OPTIMIZE | 丙烯酰胺 | （设计弹头需用`warhead_smarts`豁免，见下） |
| `diaminobenzene_1_4`/`1_2` | 氧化为醌二亚胺，遗传毒性 | REJECT | 对/邻苯二胺 | 单侧酰化的二胺 |
| `methylheterocycle_alpha_n` | 醛氧化酶(AO)代谢位点，种属差异大 | OPTIMIZE | 2-甲基吡啶 | 3-甲基吡啶(非α位) |

### 关键细节：奥希替尼不应该被这些规则误报，但确实验证到了两处真实的误报并修复

开发过程中发现 `aniline_unmasked` 和 `diaminobenzene_1_4/1_2` 最初的SMARTS
会把奥希替尼自身的结构错误标记为风险：

1. 奥希替尼的嘧啶基氨基(2-氨基嘧啶型铰链氢键基序)是激酶抑制剂里标准的、
   临床已大量验证安全的设计元素，不应该被当成"未掩蔽苯胺"——修复方式是给
   `aniline_unmasked` 加了2-氨基嘧啶/2-氨基吡啶/苯并咪唑的豁免子结构
   (`exempt_if`)，并用吉非替尼交叉验证了豁免同样生效
   (`tests/test_structural_alerts.py::test_gefitinib_aminoquinazoline_not_flagged_as_unmasked_aniline`)。
2. 奥希替尼中心苯环上有两个含氮取代基，但其中一个是酰胺(丙烯酰胺弹头上的N)，
   不满足"两个游离胺"的条件——修复方式是给两条diaminobenzene规则都加上
   `!$(N[C,S]=[O,S,N])` 排除子式(排除酰胺型氮)。

### 共价弹头豁免

路线B/C的候选分子，设计好的共价弹头本身会命中 `michael_acceptor_offtarget`
——这是预期行为，不是缺陷。调用时传入弹头SMARTS做豁免，**不要去改这条规则
本身**（这条原则直接写在yaml的advice字段里）：

```python
run_structural_alerts(mol, warhead_smarts="[CX3]=[CX3][CX3]=[OX1]")
```

### hERG风险的"子结构+计算属性"组合规则

不适合用纯SMARTS表达的规则(hERG风险、CYP3A4 TDI)单独在 `core/admet.py` 里实现
为"子结构+计算属性"组合判断：

```python
estimate_herg_risk_proxy(mol)
# 强碱性脂肪叔胺 + cLogP > 3.5 -> "high"
# 强碱性脂肪叔胺 + cLogP > 2.0 -> "medium"
# 否则 -> "low"
```

这不是hERG IC50预测模型，只是方向性的早期预警标记。

## 三态判定

`AlertReport.verdict`：任意一条REJECT命中则整体REJECT；否则任意一条
OPTIMIZE_NEEDED命中则整体OPTIMIZE_NEEDED；都没有则PASS。`workflows/run_discovery_round.py`
用这个verdict做硬门槛，REJECT的候选直接从下一步的Pareto/decision阶段排除
（但仍保留在 `rejected_pool` 里供环节8的"假阴性复活"机制使用）。
