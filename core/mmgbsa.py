"""
环节4 L3（后半）：MM-GBSA/MM-PBSA 结合自由能估算

同样需要MD轨迹(OpenMM)作为输入，本环境无法产出。这里只给出正确的接口形状和
误差量级标注(pipeline.yaml: L3误差±1.3 kcal/mol)，不猜数字。
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass
class MMGBSAResult:
    ligand_id: str
    delta_g_kcal_mol: float
    error_estimate_kcal_mol: float = 1.3  # 对应pipeline.yaml funnel.L3.error_kcal_mol
    n_frames_used: int = 0
    method: str = "mmgbsa"


def run_mmgbsa_stub(ligand_id: str, trajectory_path: str, n_frames: int = 100) -> dict:
    """占位：真实运行需要已经产出的MD轨迹(见core/md_stability.py)+ MMPBSA.py或等效工具。"""
    return {
        "ok": False,
        "reason": "需要真实MD轨迹做MM-GBSA后处理，本环境没有轨迹文件(见core/md_stability.py的说明)",
        "ligand_id": ligand_id,
    }


if __name__ == "__main__":
    print(run_mmgbsa_stub("DEMO-001", "trajectory.dcd"))
