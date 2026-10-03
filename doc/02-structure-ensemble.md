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

## 2026-10更新：突变建模不再是占位，MD平衡/构象聚类仍然是

`mutate_residue()`现在是真实函数(用PyRosetta做点突变+局部repacking，不是
"把残基名字段直接改掉"这种会产出物理上不合理结构的伪实现)，lazy import——
主.venv(Python 3.9)装不了PyRosetta，调用时诚实返回`ok=False`；真实调用需要
`.venv310/bin/python`(安装方式见`scripts/setup_docking_env.sh`或
`doc/routes/route-c-4th-gen-tki.md`里的记录)。已经在路线C上真实用过一次：
把6LUD(三突变体)的T790M突变回野生型T790，得到route C primary genotype
(L858R/C797S)的第一个真实计算结构，验证过能量(-202.5，比原结构更稳定)和
突变后身份(790=THR/797=SER/858=ARG全部确认)，详见该路线文档。

```python
from core.structures import mutate_residue

mutate_residue("6lud_receptor_only.pdb", chain="A", resnum=790, new_aa_one_letter="T",
                out_path="L858R_C797S_model.pdb", pack_radius=8.0)
# {"ok": True, "original_resname": "MET", "new_resname": "THR",
#  "total_score_after_repacking": -202.5, "out_path": "..."}
```

MD平衡/构象聚类这两步仍然是诚实占位，需要 OpenMM/GROMACS+GPU（多副本MD平衡，
典型配置是20ns×3副本，见 `config/pipeline.yaml` 的
`funnel.L3.md_length_ns`/`md_replicas`）以及真实轨迹文件（RMSD聚类）：

```python
run_md_equilibration_stub(structure_path, length_ns=500, n_replicas=3)
# -> {"ok": False, "reason": "需要OpenMM/GROMACS + GPU跑3条500ns独立轨迹..."}

cluster_representative_stub(trajectory_id, occupancy_threshold_pct=5.0)
# -> {"ok": False, "reason": "需要真实MD轨迹做RMSD聚类..."}
```

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
