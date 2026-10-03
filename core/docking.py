"""
环节4 L1-L2：快速对接 + 精细对接重打分

!! 本模块依赖的 meeko / vina 在本仓库的开发环境里装不上，见 requirements.txt 的说明 !!
  - meeko 需要 Python >= 3.10（本环境是3.9.6，meeko内部用了match语句）
  - vina 需要系统预装 Boost C++ 库
所以这里的函数都是"结构正确、可以读代码审查、但在当前环境跑不起来"的状态——
用 lazy import 包起来，import 这个模块本身不会报错，只有真正调用函数、且环境里
确实没装对应库时才会得到一个说明清楚的失败结果，不会半途崩溃也不会假装成功。

这个模块存在的第二个目的：**修正 egfr-pipline 项目里的一个真实bug**。
那个项目的 prepare_receptor_pdbqt() 实际是 `cp receptor.pdb receptor.pdbqt`——
PDB文件没有AutoDock需要的原子类型/部分电荷/可旋转键信息，只是换了个后缀名，
导致那个项目产出的全部对接结果 binding_affinity 都是 0.0（不是真实Vina打分）。
这里改用 meeko 官方推荐的 `mk_prepare_receptor.py` CLI 做真正的受体PDBQT转换。
"""
from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path


def _check_tool_available(cmd: str) -> bool:
    try:
        subprocess.run([cmd, "--help"], capture_output=True, timeout=10)
        return True
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def prepare_receptor_pdbqt(receptor_pdb: str | Path, out_pdbqt: str | Path) -> dict:
    """
    用 meeko 的 mk_prepare_receptor.py CLI 做真正的受体PDBQT转换(原子类型+部分电荷+
    可旋转键)，不是简单复制文件。要求：`pip install meeko` 且其安装的脚本在 PATH 上
    (标准 pip 安装会自动放到虚拟环境的 bin/ 目录下)。
    """
    receptor_pdb = Path(receptor_pdb)
    out_pdbqt = Path(out_pdbqt)

    if not _check_tool_available("mk_prepare_receptor.py"):
        return {
            "ok": False,
            "reason": (
                "mk_prepare_receptor.py 不在PATH上(meeko未安装或环境不兼容，见本模块顶部声明)。"
                "本环境是Python 3.9.6，meeko需要>=3.10，请用新环境重新安装: pip install meeko"
            ),
        }

    try:
        result = subprocess.run(
            [
                "mk_prepare_receptor.py",
                "--read_pdb", str(receptor_pdb),
                "-o", str(out_pdbqt.with_suffix("")),  # meeko会自己加.pdbqt后缀
                "--allow_bad_res",
            ],
            capture_output=True,
            text=True,
            timeout=300,
        )
        if result.returncode != 0:
            return {"ok": False, "reason": f"mk_prepare_receptor.py 失败: {result.stderr}"}
        return {"ok": True, "out_pdbqt": str(out_pdbqt)}
    except Exception as exc:
        return {"ok": False, "reason": str(exc)}


def prepare_ligand_pdbqt(ligand_pdb: str | Path, out_pdbqt: str | Path) -> dict:
    """用 meeko.MoleculePreparation 把单个配体的3D结构(.pdb)转成Vina需要的.pdbqt。"""
    try:
        from meeko import MoleculePreparation
    except ImportError as exc:
        return {"ok": False, "reason": f"meeko 未安装或不兼容当前Python版本: {exc}"}

    from rdkit import Chem

    mol = Chem.MolFromPDBFile(str(ligand_pdb), removeHs=False)
    if mol is None:
        return {"ok": False, "reason": f"RDKit无法解析配体结构: {ligand_pdb}"}

    preparator = MoleculePreparation()
    preparator.prepare(mol)
    pdbqt_string = preparator.write_pdbqt_string()
    Path(out_pdbqt).write_text(pdbqt_string)
    return {"ok": True, "out_pdbqt": str(out_pdbqt)}


@dataclass
class DockingResult:
    ligand_id: str
    best_affinity_kcal_mol: float | None
    n_poses: int
    docked_pdbqt_path: str | None
    ok: bool
    error: str | None = None


def run_vina_docking(
    receptor_pdbqt: str | Path,
    ligand_pdbqt: str | Path,
    center: tuple[float, float, float],
    box_size: tuple[float, float, float],
    out_pdbqt: str | Path,
    exhaustiveness: int = 8,
    n_poses: int = 3,
) -> DockingResult:
    """
    环节4 L1: AutoDock Vina 快速对接。center/box_size 来自 core/structures.py 准备的
    结构配合已知共晶配体定义的口袋坐标(不要凭空猜)。
    结果只用于排除(pipeline.yaml 的 l1_exclude_above_kcal_mol)，不要直接拿Vina分数排名——
    见设计文档环节4开头"用Vina做最终排序，等于用卷尺量头发丝"。
    """
    try:
        from vina import Vina
    except ImportError as exc:
        return DockingResult(
            ligand_id=str(ligand_pdbqt),
            best_affinity_kcal_mol=None,
            n_poses=0,
            docked_pdbqt_path=None,
            ok=False,
            error=f"vina 未安装或不兼容当前Python版本: {exc}（还需要系统预装Boost库）",
        )

    try:
        v = Vina(sf_name="vina")
        v.set_receptor(str(receptor_pdbqt))
        v.set_ligand_from_file(str(ligand_pdbqt))
        v.compute_vina_maps(center=list(center), box_size=list(box_size))
        v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
        energies = v.energies()
        best_affinity = float(energies[0][0])
        v.write_poses(str(out_pdbqt), n_poses=1, overwrite=True)
        return DockingResult(
            ligand_id=str(ligand_pdbqt),
            best_affinity_kcal_mol=round(best_affinity, 2),
            n_poses=n_poses,
            docked_pdbqt_path=str(out_pdbqt),
            ok=True,
        )
    except Exception as exc:
        return DockingResult(
            ligand_id=str(ligand_pdbqt), best_affinity_kcal_mol=None, n_poses=0, docked_pdbqt_path=None, ok=False, error=str(exc)
        )


if __name__ == "__main__":
    print("=== 受体准备 (预期: ok=False，本环境meeko不可用) ===")
    print(prepare_receptor_pdbqt("receptor_clean.pdb", "receptor_clean.pdbqt"))

    print("\n=== Vina对接 (预期: ok=False，本环境vina不可用) ===")
    r = run_vina_docking(
        receptor_pdbqt="receptor_clean.pdbqt",
        ligand_pdbqt="ligand.pdbqt",
        center=(-14.2, 33.5, 22.8),
        box_size=(20.0, 20.0, 20.0),
        out_pdbqt="docked.pdbqt",
    )
    print(r)
    print(
        "\n两个函数都诚实返回了 ok=False + 具体原因，而不是像egfr-pipline那样"
        "静默产出一个看起来成功、实际是0.0 kcal/mol的假结果。"
    )
