# 环节3：类似物枚举

代码：[`core/enumerate.py`](../core/enumerate.py) · 配置：[`config/rgroup_libraries/`](../config/rgroup_libraries/) · 测试：[`tests/test_enumerate.py`](../tests/test_enumerate.py) [`tests/test_rgroup_libraries.py`](../tests/test_rgroup_libraries.py) · **状态：真实可跑**

## 设计思路：按口袋区域定向枚举，不是随机加官能团

原设计文档"环节3.2"强调：有效的R基团枚举不是对着一个骨架随机尝试官能团，而是
按每个可变位点在蛋白口袋里对应的具体区域（铰链区/797入口区/疏水背袋/溶剂暴露区）
去挑选在**那个区域**有物理意义的片段。比如溶剂暴露区决定PK和hERG风险，所以
`config/rgroup_libraries/pocket_regions.yaml` 里这个区域刻意富集了弱碱性/中性
的替代基团（吗啉、N-甲基哌嗪、砜化哌嗪、氟代氮杂环丁烷等），避免"强碱性叔胺+
高亲脂性"这个hERG风险的经典组合。

## 数据结构

**骨架模板**（`config/rgroup_libraries/scaffolds.yaml`）：用同位素标记的哑原子
`[1*]`/`[2*]`/`[3*]`标记可变位点，每个位点声明对应哪个口袋区域：

```yaml
templates:
  - id: "demo_aminopyrimidine_biphenyl"
    status: "illustrative_only"   # 示例性骨架，不是真实候选起始点
    core_smiles: "[1*]c1ccc(Nc2nccc(-c3ccc([2*])cc3[3*])n2)cc1"
    attachment_points:
      - {label: R1, isotope: 1, region: entrance_797}
      - {label: R2, isotope: 2, region: hydrophobic_back_pocket}
      - {label: R3, isotope: 3, region: solvent_exposed}
```

**片段库**（`config/rgroup_libraries/pocket_regions.yaml`）：四个口袋区域
（hinge/entrance_797/hydrophobic_back_pocket/solvent_exposed），每个区域下是
一组用裸 `*` 标记连接点的SMILES片段。

## 拼接机制：RDKit molzip + 同位素标记

片段库里的 `*` 统一不带同位素标记（同一份库要能复用到不同骨架的不同位点上，
不能在配置里写死某个编号）。拼接时需要先把片段的哑原子同位素改成骨架要求的
编号，再用 `Chem.molzip` 按同位素配对连接：

```python
def _mol_with_isotope_dummy(smiles, force_isotope=None):
    mol = Chem.MolFromSmiles(smiles)
    if force_isotope is not None:
        dummy_atoms = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
        dummy_atoms[0].SetIsotope(force_isotope)
    return mol
```

## 开发时真实踩过的一个坑（已修复，已加回归测试）

**如果不做上面这步同位素改写，直接把裸 `*` 片段和带编号的骨架喂给
`molzip`，RDKit既不报错也不抛异常，而是悄悄把片段原样留成游离的
`.`分子**——产出类似 `*C.*N1CCOCC1.[1*]c1ccc(...)` 这样的多组分SMILES，
`Chem.SanitizeMol()` 对这种结构完全不会报错，表面上看起来是一个合法的分子。
这个bug是在检查实际生成的SMILES字符串（注意到有可疑的`.`分隔符）时才发现的，
而不是靠"代码跑起来没报错"发现的。

修复方式：
1. `assemble()` 里强制用 `force_isotope` 把片段同位素改成骨架要求的编号；
2. 拼接后做显式连通性校验：`len(Chem.GetMolFrags(result)) != 1` 就返回
   `None`，宁可整条候选丢弃也不能让断裂分子混进枚举结果；
3. 额外检查拼接后是否还有残留哑原子（说明某个连接点没被消耗）。

`tests/test_enumerate.py::test_assemble_rejects_mismatched_isotope_as_broken_not_silent`
专门复现了这个坏路径，确保"molzip不报错≠拼接成功"这个事实被测试锁定，不会
在未来的重构里被意外引入。

## 接口

```python
from core.enumerate import enumerate_from_scaffold

df = enumerate_from_scaffold("demo_aminopyrimidine_biphenyl")
# df.columns -> ['smiles', 'scaffold_id', 'R1_name', 'R2_name', 'R3_name']
```

demo骨架的三个位点分别对应entrance_797(6个片段)、hydrophobic_back_pocket(6个)、
solvent_exposed(9个)，全组合枚举产出 6×6×9 = 324 个连通、去重后的候选分子
（已验证，见 `tests/test_enumerate.py::test_enumerate_count_matches_cartesian_product_of_region_sizes`）。

真实商业R基团库通常是成百上千个/区域，全组合会爆炸，这种规模要用
`max_combinations` 参数做随机抽样封顶，或者换成原设计文档"环节3.3"描述的
"受约束生成"而不是穷举——这部分尚未实现，当前 `enumerate_from_scaffold()`
只支持全组合枚举。
