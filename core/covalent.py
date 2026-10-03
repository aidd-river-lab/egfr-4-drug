"""
环节4.3：共价对接（仅chemistry_route=covalent_new_site/covalent_pan_mutant_broad时启用，
见 target_profile.yaml 的 chemistry_route.covalent_warhead_enabled）

真实共价对接需要 Schrödinger CovDock / AutoDock4-covalent / Rosetta，本环境未安装。
真实、可测试的部分：Bürgi-Dunitz攻击角判据——这是纯几何判断，不需要真的跑过
共价对接就能验证逻辑对不对（设计文档环节4.3第2步："攻击角约105°"）。

共价药的核心认知（环节4.3末尾）：**IC50没有意义，要看kinact/KI**。KI(可逆亲和力)和
kinact(成键速率)优化方向不同——提高KI靠结合位姿，提高kinact靠弹头几何对齐和
亲核体pKa。所以 CovalentDockingResult 把两者分开存，不合并成一个数。
"""
from __future__ import annotations

from dataclasses import dataclass

import yaml
from pathlib import Path

from core.route_config import route_config_dir

BD_ANGLE_TOLERANCE_DEGREE = 15  # 允许的偏差范围，105°±15°


@dataclass
class CovalentDockingResult:
    ligand_id: str
    nucleophile_residue: str  # 比如 "Lys745" / "Cys797"(如果走的是经典弹头路线)
    warhead_carbon_to_nucleophile_distance_angstrom: float
    attack_angle_degree: float
    ki_nm: float | None = None  # 可逆亲和力，决定选择性(环节5.2)
    kinact_per_second: float | None = None  # 成键速率，决定共价效率
    gsh_half_life_hours: float | None = None  # 弹头本征反应性

    @property
    def kinact_over_ki(self) -> float | None:
        """共价药的主终点：kinact/KI，不是IC50。"""
        if self.ki_nm is None or self.kinact_per_second is None or self.ki_nm == 0:
            return None
        return self.kinact_per_second / self.ki_nm


def evaluate_attack_geometry(result: CovalentDockingResult, config_dir: Path | None = None) -> dict:
    """
    判断预结合位姿的弹头几何是否合理：攻击角接近Bürgi-Dunitz角(~105°)。
    这条判据是纯几何检查，用合成数据就能完整验证（见 __main__）。

    config_dir: 不传则默认路线C的pipeline.yaml；路线A/B传各自的 route_config_dir(...)。
    """
    config_dir = config_dir or route_config_dir()
    with open(config_dir / "pipeline.yaml", encoding="utf-8") as f:
        target_angle = yaml.safe_load(f)["covalent_route"]["attack_angle_bd_degree"]

    deviation = abs(result.attack_angle_degree - target_angle)
    geometry_ok = deviation <= BD_ANGLE_TOLERANCE_DEGREE
    return {
        "ligand_id": result.ligand_id,
        "target_angle_degree": target_angle,
        "actual_angle_degree": result.attack_angle_degree,
        "deviation_degree": round(deviation, 1),
        "geometry_ok": geometry_ok,
        "kinact_over_ki": result.kinact_over_ki,
        "note": "共价对接的IC50没有意义；排序/决策请用kinact_over_ki，KI单独拿出来做选择性分析(环节4.3末尾)",
    }


def run_covalent_docking_stub(ligand_id: str, nucleophile_residue: str) -> dict:
    """占位：真实运行需要Schrödinger CovDock/AutoDock4-covalent/Rosetta，本环境未安装。"""
    return {
        "ok": False,
        "reason": "需要CovDock/AutoDock4-covalent/Rosetta做共价对接，本环境未安装，不编造攻击角/距离数值",
        "ligand_id": ligand_id,
        "nucleophile_residue": nucleophile_residue,
    }


if __name__ == "__main__":
    print("=== run_covalent_docking_stub (预期 ok=False) ===")
    print(run_covalent_docking_stub("DEMO-001", "Lys745"))

    print("\n=== evaluate_attack_geometry 逻辑验证(合成数据) ===")
    good = CovalentDockingResult(
        ligand_id="DEMO-GOOD", nucleophile_residue="Lys745",
        warhead_carbon_to_nucleophile_distance_angstrom=3.2, attack_angle_degree=103,
        ki_nm=50, kinact_per_second=0.02,
    )
    bad = CovalentDockingResult(
        ligand_id="DEMO-BAD", nucleophile_residue="Lys745",
        warhead_carbon_to_nucleophile_distance_angstrom=5.1, attack_angle_degree=60,
        ki_nm=800, kinact_per_second=0.001,
    )
    print(evaluate_attack_geometry(good))
    print(evaluate_attack_geometry(bad))
