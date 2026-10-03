"""
环节2：结构生态系统（结构库 structure ensemble 管理）

真正能在没有GPU/MD软件的开发环境里跑通的部分：
  - fetch_pdb() / clean_chain()：从RCSB下载 + 按链清洗，这是环节2.2第1步，纯IO操作
这个模块里"需要MD/突变建模软件"的部分（环节2.2第2-4步）用明确标注"未执行"的
orchestration接口表示，不假装算出了构象聚类结果——原因和 core/synthesis.py 的
run_retrosynthesis_stub 一样：没有真实计算后端时，诚实地说"没跑"比编一个假结果更有用。
"""
from __future__ import annotations

import urllib.request
from dataclasses import dataclass, field
from pathlib import Path


def fetch_pdb(pdb_id: str, out_dir: str | Path) -> Path:
    """从 RCSB PDB 下载结构，返回本地文件路径。真实网络IO，会实际落盘。"""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{pdb_id.lower()}_raw.pdb"
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    urllib.request.urlretrieve(url, out_path)
    return out_path


def clean_chain(raw_pdb_path: str | Path, chain: str, out_path: str | Path, keep_hetatm: bool = False) -> dict:
    """
    从原始PDB里只保留指定链的标准蛋白原子(ATOM记录)。
    keep_hetatm=True 时额外保留该链的HETATM记录(配体/辅因子/结构水)——
    egfr-pipline 那版实现默认丢弃了所有HETATM，如果后续要做"保留共晶配体用于口袋定义"
    或"保留结构水"(环节2.1里明确提到的关键水分子)，需要把这个开关打开。
    """
    raw_pdb_path = Path(raw_pdb_path)
    out_path = Path(out_path)
    n_atom = 0
    n_hetatm = 0
    with open(raw_pdb_path) as f_in, open(out_path, "w") as f_out:
        for line in f_in:
            record_type = line[:6].strip()
            if len(line) < 22:
                continue
            chain_col = line[21]
            if record_type == "ATOM" and chain_col == chain:
                f_out.write(line)
                n_atom += 1
            elif keep_hetatm and record_type == "HETATM" and chain_col == chain:
                f_out.write(line)
                n_hetatm += 1
    return {"out_path": str(out_path), "n_atom_lines": n_atom, "n_hetatm_lines": n_hetatm}


@dataclass
class StructureEnsembleMember:
    """环节2.4 要求的结构库元数据字段——每个构象都要能回答"这是从哪来的、怎么做出来的"。"""

    genotype: str  # 对应 target_profile.yaml 里的基因型id，比如 del19_C797S / EGFR_WT
    source_pdb: str
    mutation_method: str | None = None  # 比如 "Rosetta residue mutation" / "none(直接来自实验结构)"
    md_trajectory_id: str | None = None
    cluster_occupancy_pct: float | None = None
    prep_date: str | None = None
    forcefield_version: str | None = None
    local_path: str | None = None
    notes: str = ""


def mutate_residue(
    structure_path: str, chain: str, resnum: int, new_aa_one_letter: str, out_path: str, pack_radius: float = 8.0
) -> dict:
    """
    环节2.2第2步：真实的点突变建模(2026-10起不再是纯占位)。

    用 PyRosetta 的 mutate_residue()：换残基身份 + 以pack_radius(Å)为半径对周围
    侧链做局部repacking，不是简单的"把PDB文本里的残基名字段改掉"那种会产出物理上
    不合理结构(侧链和周围原子重叠)的伪实现。

    主.venv(Python 3.9)装不了PyRosetta，这里用lazy import——import这个模块本身
    不会报错，只有真正调用且当前解释器没装pyrosetta时才诚实返回ok=False。
    真实调用必须用`.venv310/bin/python`(见 scripts/setup_docking_env.sh，
    PyRosetta的安装方式见该脚本或doc/routes/route-c-4th-gen-tki.md的记录——
    注意pyrosetta-installer内部调用裸`pip`而不是`sys.executable -m pip`，
    如果PATH上别的Python版本的pip排在前面会装错环境，需要显式把目标venv的
    bin目录放到PATH最前面)。

    new_aa_one_letter: 单字母氨基酸代码(PyRosetta的约定，不是三字母)。
    """
    try:
        import pyrosetta
        from pyrosetta.toolbox.mutants import mutate_residue as _pyrosetta_mutate
    except ImportError as exc:
        return {
            "ok": False,
            "reason": f"需要PyRosetta做突变建模+局部repacking，当前解释器未安装: {exc}",
            "structure_path": structure_path,
            "resnum": resnum,
            "new_aa_one_letter": new_aa_one_letter,
        }

    pyrosetta.init("-mute all")
    pose = pyrosetta.pose_from_pdb(structure_path)
    pdb_info = pose.pdb_info()
    pose_idx = pdb_info.pdb2pose(chain, resnum)
    if pose_idx == 0:
        return {"ok": False, "reason": f"链{chain}里找不到PDB残基编号{resnum}", "structure_path": structure_path}

    original_resname = pose.residue(pose_idx).name3()
    _pyrosetta_mutate(pose, pose_idx, new_aa_one_letter, pack_radius=pack_radius)
    new_resname = pose.residue(pose_idx).name3()

    scorefxn = pyrosetta.get_fa_scorefxn()
    score_after = scorefxn(pose)

    pose.dump_pdb(out_path)
    return {
        "ok": True,
        "chain": chain,
        "resnum": resnum,
        "original_resname": original_resname,
        "new_resname": new_resname,
        "total_score_after_repacking": round(score_after, 1),
        "out_path": out_path,
    }


def run_md_equilibration_stub(structure_path: str, length_ns: int, n_replicas: int) -> dict:
    """
    环节2.2第3步：多副本MD平衡的占位接口，真实需要 OpenMM/GROMACS + GPU，耗时数小时到数天
    （见 doc/02-structure-ensemble.md 的算力预算）。本环境没有GPU，不编造轨迹文件。
    """
    return {
        "ok": False,
        "reason": f"需要OpenMM/GROMACS + GPU跑{n_replicas}条{length_ns}ns独立轨迹，本环境无GPU，未执行",
        "structure_path": structure_path,
    }


def cluster_representative_stub(trajectory_id: str, occupancy_threshold_pct: float = 5.0) -> dict:
    """环节2.2第4步：RMSD聚类/马尔可夫态模型，需要真实轨迹文件，本环境没有可用轨迹。"""
    return {
        "ok": False,
        "reason": "需要真实MD轨迹做RMSD聚类，见run_md_equilibration_stub，本环境未产出轨迹",
        "trajectory_id": trajectory_id,
    }


if __name__ == "__main__":
    import tempfile

    with tempfile.TemporaryDirectory() as tmpdir:
        print("=== 下载 6S9C (EGFR C797S 共晶结构) ===")
        raw_path = fetch_pdb("6S9C", tmpdir)
        print(f"下载完成: {raw_path}, 大小 {raw_path.stat().st_size} bytes")

        print("\n=== 清洗A链 ===")
        clean_path = Path(tmpdir) / "receptor_A_clean.pdb"
        result = clean_chain(raw_path, chain="A", out_path=clean_path)
        print(result)

        member = StructureEnsembleMember(
            genotype="del19_C797S",
            source_pdb="6S9C",
            mutation_method="none(直接来自实验共晶结构)",
            local_path=str(clean_path),
            notes="演示用途；真实项目需要再跑MD平衡+聚类才能进对接(见环节2.2第3-4步)",
        )
        print(f"\n=== 结构库元数据记录 ===\n{member}")

        print("\n=== MD/突变占位接口（预期诚实返回 ok=False）===")
        print(run_md_equilibration_stub(str(clean_path), length_ns=500, n_replicas=3))
