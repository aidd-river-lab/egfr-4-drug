"""
环节4.3：共价对接（仅chemistry_route=covalent_new_site/covalent_pan_mutant_broad时启用，
见 target_profile.yaml 的 chemistry_route.covalent_warhead_enabled）

真实共价对接需要 Schrödinger CovDock / AutoDock4-covalent / Rosetta，本环境未安装。
真实、可测试的部分：Bürgi-Dunitz攻击角判据——这是纯几何判断，不需要真的跑过
共价对接就能验证逻辑对不对（设计文档环节4.3第2步："攻击角约105°"）。

共价药的核心认知（环节4.3末尾）：**IC50没有意义，要看kinact/KI**。KI(可逆亲和力)和
kinact(成键速率)优化方向不同——提高KI靠结合位姿，提高kinact靠弹头几何对齐和
亲核体pKa。所以 CovalentDockingResult 把两者分开存，不合并成一个数。

2026-10补充：路线C真实候选(demo_aminopyrimidine_biphenyl骨架)实际上是
reversible_orthosteric(chemistry_route.primary，covalent_warhead_enabled=false)，
不含共价弹头——evaluate_attack_geometry()这个判据函数从未在真实候选上跑过，
因为真实候选根本没有弹头可供测量。真正有弹头的只有奥希替尼这个外部参照分子。
`extract_michael_acceptor_geometry()`补上了"从一个真实对接姿态里量出弹头反应碳
到亲核原子的距离/角度"这一步(之前只有手填数值进`CovalentDockingResult`的演示，
没有从真实3D结构提取过)，用奥希替尼+真实PDB 4ZAU(AZD9291与野生型EGFR复合物，
Cys797完整保留)验证过一次，见`doc/15`和`scripts/run_covalent_geometry_check.py`。
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import yaml

from core.route_config import route_config_dir

BD_ANGLE_TOLERANCE_DEGREE = 15  # 允许的偏差范围，105°±15°

_PDBQT_ELEMENT_MAP = {
    "OA": "O", "NA": "N", "SA": "S", "HD": "H", "A": "C", "N": "N",
    "C": "C", "O": "O", "S": "S", "H": "H", "F": "F", "Cl": "Cl",
    "Br": "Br", "I": "I", "P": "P",
}


def _parse_pdbqt_atoms(pdbqt_path: str | Path) -> list[dict]:
    """解析PDBQT的ATOM/HETATM行：标准PDB列(resname/resnum/坐标) + 末尾AutoDock原子类型。
    不依赖CONECT/成键信息(PDBQT里这部分不可靠)，只拿坐标和元素，连通性后面用纯几何推断。"""
    atoms = []
    with open(pdbqt_path) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            autodock_type = line[77:79].strip() or line[76:78].strip()
            element = _PDBQT_ELEMENT_MAP.get(autodock_type, autodock_type[:1])
            atoms.append(
                {
                    "atom_name": line[12:16].strip(),
                    "resname": line[17:20].strip(),
                    "resnum": int(line[22:26].strip()),
                    "chain": line[21],
                    "x": float(line[30:38]),
                    "y": float(line[38:46]),
                    "z": float(line[46:54]),
                    "element": element,
                }
            )
    return atoms


def _find_unique_michael_beta_carbon(ligand_atoms: list[dict]) -> dict | None:
    """纯几何规则识别丙烯酰胺端位烯烃碳(Cβ)：degree-1重原子邻居是碳、键长<1.40Å
    (C=C双键长度量级，区别于C-N单键~1.45-1.47Å)。已用RDKit在奥希替尼真实3D结构上
    验证过这个组合在全分子范围内是唯一的(core/covalent.py模块文档记录的验证过程，
    见scripts/run_covalent_geometry_check.py)。找不到恰好1个符合条件的碳，诚实返回
    None，不强行凑一个。"""
    heavy = [a for a in ligand_atoms if a["element"] != "H"]
    candidates = []
    for atom in heavy:
        if atom["element"] != "C":
            continue
        neighbors = []
        for other in heavy:
            if other is atom:
                continue
            dist = math.dist((atom["x"], atom["y"], atom["z"]), (other["x"], other["y"], other["z"]))
            if dist < 1.70:  # 重原子间成键距离的宽松上限
                neighbors.append((other, dist))
        if len(neighbors) == 1 and neighbors[0][0]["element"] == "C" and neighbors[0][1] < 1.40:
            candidates.append({"beta": atom, "alpha": neighbors[0][0], "bond_length": neighbors[0][1]})
    if len(candidates) != 1:
        return None
    return candidates[0]


def extract_michael_acceptor_geometry(
    docked_ligand_pdbqt: str | Path,
    receptor_pdbqt: str | Path,
    nucleophile_resname: str,
    nucleophile_resnum: int,
    nucleophile_atom_name: str,
) -> dict:
    """
    从一个真实对接姿态(非共价，Vina给不出真正的共价键，只能给反应前的encounter
    complex姿态)里，量出丙烯酰胺弹头端位烯烃碳(Cβ)到受体亲核原子(比如Cys-SG)的
    真实距离，以及攻击角。

    攻击角定义(本仓库自己采用的工作定义，不是某个标准共价对接软件如Schrödinger
    CovDock的精确复刻，`core/covalent.py`本身没有预先规定用哪几个原子算这个角，
    这里显式写清楚)：angle(亲核原子, Cβ, Cα)，Cα是Cβ唯一的重原子邻居(双键另一端)。
    如果要对齐某个具体共价对接软件的精确角度惯例，需要药物化学/结构生物学专家确认。

    找不到唯一的丙烯酰胺端位碳(比如这个配体根本不含迈克尔受体弹头)，或者受体里
    找不到指定的亲核原子，诚实返回ok=False，不编造数值。
    """
    ligand_atoms = _parse_pdbqt_atoms(docked_ligand_pdbqt)
    beta_alpha = _find_unique_michael_beta_carbon(ligand_atoms)
    if beta_alpha is None:
        return {
            "ok": False,
            "reason": "配体里没有找到唯一的丙烯酰胺端位烯烃碳(Cβ)——这个分子可能不含"
            "迈克尔受体弹头，或者含有多个类似基团导致无法唯一确定",
        }

    receptor_atoms = _parse_pdbqt_atoms(receptor_pdbqt)
    nucleophile_matches = [
        a
        for a in receptor_atoms
        if a["resname"] == nucleophile_resname
        and a["resnum"] == nucleophile_resnum
        and a["atom_name"] == nucleophile_atom_name
    ]
    if len(nucleophile_matches) != 1:
        return {
            "ok": False,
            "reason": f"受体里找到{len(nucleophile_matches)}个{nucleophile_resname}"
            f"{nucleophile_resnum}-{nucleophile_atom_name}原子，需要恰好1个",
        }
    nucleophile = nucleophile_matches[0]

    beta, alpha = beta_alpha["beta"], beta_alpha["alpha"]
    p_nu = np.array([nucleophile["x"], nucleophile["y"], nucleophile["z"]])
    p_beta = np.array([beta["x"], beta["y"], beta["z"]])
    p_alpha = np.array([alpha["x"], alpha["y"], alpha["z"]])

    distance = float(np.linalg.norm(p_nu - p_beta))
    v1 = p_nu - p_beta
    v2 = p_alpha - p_beta
    cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    angle_degree = float(np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0))))

    return {
        "ok": True,
        "nucleophile": f"{nucleophile_resname}{nucleophile_resnum}-{nucleophile_atom_name}",
        "beta_carbon_to_nucleophile_distance_angstrom": round(distance, 2),
        "attack_angle_degree": round(angle_degree, 1),
        "beta_alpha_bond_length_angstrom": round(beta_alpha["bond_length"], 3),
        "note": "这是非共价对接给出的反应前encounter complex几何，不是共价键形成后的键长"
        "(真实共价键长约1.8Å)——Vina不能模拟反应本身，只能预测反应前弹头有没有摆在"
        "合理的攻击轨迹上",
    }


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
