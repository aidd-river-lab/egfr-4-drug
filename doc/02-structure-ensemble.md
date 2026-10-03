# 环节2：结构系综准备

代码：[`core/structures.py`](../core/structures.py) · **状态：下载/清洗部分真实可跑，MD/突变建模部分诚实占位**

## 为什么需要"系综"而不是单一晶体结构

对接/MM-GBSA/FEP都是对着一个固定的蛋白构象算结合能，但蛋白是动态的——单一
晶体结构的侧链构象、水分子位置、口袋开合状态都只是无数个可能状态里的一个快照。
如果只用一个构象做所有后续打分，很容易出现"分子A在这个快照里对接分数很差，
但实际上它能稳定结合另一个（没被采样到的）构象"这种假阴性。

环节2的目标是为每个基因型准备一组有代表性的构象（"系综"），而不是单一静态结构。

## 真实可跑的部分：下载 + 清洗

```python
from core.structures import fetch_pdb, clean_chain

raw_path = fetch_pdb("6S9C", "/tmp/structures")       # 从RCSB真实下载
result = clean_chain(raw_path, chain="A", out_path="receptor_A.pdb")
# result = {"out_path": ..., "n_atom_lines": 2433, "n_hetatm_lines": 0}
```

`fetch_pdb()` 是真实的网络IO，直接从 `https://files.rcsb.org/download/` 下载。
`clean_chain()` 按链过滤 `ATOM`/`HETATM` 记录，纯文本处理，不依赖任何外部工具。

**与姊妹项目 `egfr-pipline` 的一个重要差异**：那个项目的清洗函数默认丢弃所有
`HETATM` 记录（配体、辅因子、结构水全部丢失）。本模块加了 `keep_hetatm` 开关，
如果后续要做"保留共晶配体定义对接口袋坐标"或"保留关键结构水"（原设计文档
环节2.1明确提到结构水对口袋定义的重要性），把这个开关打开即可。

## 诚实占位的部分：突变建模 / MD平衡 / 构象聚类

这三步需要 Rosetta/Maestro（突变建模+局部能量最小化）或 OpenMM/GROMACS+GPU
（多副本MD平衡，典型配置是20ns×3副本，见 `config/pipeline.yaml` 的
`funnel.L3.md_length_ns`/`md_replicas`）以及真实轨迹文件（RMSD聚类）。
本环境都不具备，对应函数诚实返回 `ok=False`：

```python
mutate_residue_stub(structure_path, resnum=797, new_resname="SER")
# -> {"ok": False, "reason": "需要 Rosetta/Maestro 做突变建模..."}

run_md_equilibration_stub(structure_path, length_ns=500, n_replicas=3)
# -> {"ok": False, "reason": "需要OpenMM/GROMACS + GPU跑3条500ns独立轨迹..."}

cluster_representative_stub(trajectory_id, occupancy_threshold_pct=5.0)
# -> {"ok": False, "reason": "需要真实MD轨迹做RMSD聚类..."}
```

这里**没有**做"把残基名字段直接改掉"这种文本替换式的伪突变——那样产出的结构
在物理上是不合理的（侧链会和周围原子碰撞，键长键角不对），比明确说"没做"更危险，
因为它看起来像是跑完了一步。

## 数据模型

`StructureEnsembleMember` dataclass（字段照抄原设计文档"环节2.4"的结构库元数据
表）：

```python
@dataclass
class StructureEnsembleMember:
    genotype: str                          # 如 del19_C797S / EGFR_WT
    source_pdb: str
    mutation_method: str | None = None     # "Rosetta residue mutation" / "none(实验结构)"
    md_trajectory_id: str | None = None
    cluster_occupancy_pct: float | None = None
    prep_date: str | None = None
    forcefield_version: str | None = None
    local_path: str | None = None
    notes: str = ""
```

对应落库表：`db/schema.sql` 的 `structure_ensemble`。

## 算力预算参考（来自原设计文档，接入真实MD时按此规划）

单个基因型：20ns × 3副本的MD平衡，GPU上是数小时到一天量级；本项目需要覆盖
4个主要/次要基因型 + EGFR_WT + 7个脱靶激酶，如果每个都要做系综，总算力需求是
数十到上百GPU小时。这是环节2需要排期和申请GPU资源的主要原因，也是本仓库选择
"先把接口定对，等有GPU资源再接入"而不是"在CPU上硬跑一个阉割版MD"的原因——
阉割版MD(比如只跑1ns)产出的构象多样性不具代表性，反而会给后续打分引入系统性
偏差。
