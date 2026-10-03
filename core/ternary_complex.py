"""
路线B核心难点：warhead-EGFR-E3连接酶三元复合物预测

二元对接(warhead结合EGFR、E3配体结合CRBN)只是必要条件，不是充分条件——
真正决定降解效率的是三元复合物能不能形成，以及形成时warhead/E3配体的结合界面
会不会互相干扰(负协同)还是互相促进(正协同)。这是二元对接完全捕捉不到的。

2026-10更新：不再是纯占位。放弃完整复刻PRosettaC发表论文的确切流程(需要
PatchDock+完整Rosetta C++套件+PyMOL+HPC调度系统，详见
doc/routes/route-b-degrader.md的记录)，改用PyRosetta自带的刚体对接协议
(DockingProtocol家族，不需要PatchDock/调度系统)搭一个更轻量的几何筛选版本——
原理一样(蛋白-蛋白对接 + linker几何兼容性过滤)，但不是发表论文验证过的确切方法，
准确度没有benchmark数据支撑。

真实跑过一次(路线B的warhead vs 真实CRBN结构4TZ4)：500次独立随机刚体对接+
滑入接触，exit vector距离最小只有26.8Å，而最长的linker(peg3)实测span也只有
11.9Å——**0/500次找到几何兼容的姿态**。这是一个真实、有统计量支撑的发现，不是
bug：均匀随机的刚体朝向采样，两个蛋白表面上各自很窄的exit vector朝向恰好对上
的概率天然很低，这正是真实PRosettaC论文要用PatchDock(系统性穷举表面补丁而不是
均匀随机)而不是随机采样的原因。详见doc/routes/route-b-degrader.md。

本模块里还有两类真实、可测试、不需要PyRosetta就能验证的部分：
  1. hook effect(钩状效应)判定——双功能降解剂公认的真实药理学现象，纯浓度/
     平衡常数数学关系。
  2. measure_linker_span()——纯RDKit构象生成，不需要PyRosetta。
  3. kabsch_rigid_transform()——纯numpy，把对接前后的刚体变换用在"固定坐标系
     里追踪一个不属于被对接蛋白链本身的外部参考点(比如配体的exit vector原子)"
     这个问题上，是上面500次对接实验分析用的核心数学工具。
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np


@dataclass
class TernaryComplexResult:
    """三元复合物评估结果的字段schema。真实3D预测部分留空(None)，不编造数值。"""

    warhead_id: str
    e3_ligand_id: str
    linker_id: str
    warhead_kd_nm: float | None = None  # warhead对EGFR的二元结合常数(可以从L1对接/FEP拿到)
    e3_ligand_kd_nm: float | None = None  # E3配体对CRBN的二元结合常数(通常是已知文献值)
    cooperativity_alpha: float | None = None  # 协同因子：>1正协同，<1负协同，=1无协同；需要真实三元复合物数据才能测
    predicted_dc50_nm: float | None = None
    predicted_dmax_pct: float | None = None


def evaluate_hook_effect_risk(
    warhead_kd_nm: float, e3_ligand_kd_nm: float, dose_range_nm: tuple[float, float]
) -> dict:
    """
    Hook effect判定：当给药浓度上限超过两个二元结合常数中较大者的若干倍时，
    两个二元结合会分别趋于饱和，三元复合物浓度反而下降(钟形剂量-响应曲线)。
    这里用一个保守的经验阈值(10倍)标记进入风险区间的起点，不是精确预测IC50会跌多少——
    精确预测需要完整的三元复合物平衡方程求解(还需要未知的cooperativity_alpha)，
    这里只做方向性的风险标记，诚实反映计算能力的边界。

    warhead_kd_nm/e3_ligand_kd_nm: 两个二元结合的平衡解离常数
    dose_range_nm: 计划考察的给药/实测浓度范围 (low, high)
    """
    limiting_kd = max(warhead_kd_nm, e3_ligand_kd_nm)
    hook_risk_threshold_nm = 10 * limiting_kd
    dose_low, dose_high = dose_range_nm

    in_risk_zone = dose_high >= hook_risk_threshold_nm
    safe_margin_fold = hook_risk_threshold_nm / dose_low if dose_low > 0 else None

    return {
        "limiting_kd_nm": limiting_kd,
        "hook_risk_threshold_nm": hook_risk_threshold_nm,
        "dose_range_nm": dose_range_nm,
        "hook_effect_risk": in_risk_zone,
        "safe_margin_fold": round(safe_margin_fold, 1) if safe_margin_fold else None,
        "note": (
            "剂量上限超过两个二元结合常数中较大者的10倍，进入hook effect风险区间"
            if in_risk_zone
            else "剂量范围在风险阈值以内，但这只是方向性判断，不替代真实细胞DC50/Dmax曲线"
        ),
    }


def run_ternary_complex_stub(warhead_id: str, e3_ligand_id: str, linker_id: str) -> dict:
    """
    占位：cooperativity_alpha/DC50/Dmax这几个量需要对三元复合物做能量分解(三元复合物能量
    对比两个二元复合物各自的能量)，比下面geometry screening这一步贵得多，也还没做。
    不编造这几个数值。几何筛选(有没有姿态能让linker搭上)见
    run_protein_protein_docking_trials() + evaluate_ternary_complex_geometry()。
    """
    return {
        "ok": False,
        "reason": "cooperativity_alpha/DC50/Dmax需要三元复合物能量分解，比几何筛选更贵，还没做，不编造数值",
        "warhead_id": warhead_id,
        "e3_ligand_id": e3_ligand_id,
        "linker_id": linker_id,
    }


# ------------------------------------------------------------
# 纯数学/纯RDKit部分：不需要PyRosetta，可以独立测试
# ------------------------------------------------------------
def kabsch_rigid_transform(before: np.ndarray, after: np.ndarray) -> Callable[[np.ndarray], np.ndarray]:
    """
    Kabsch算法：给定同一组点在"刚体移动前"和"刚体移动后"的坐标，求解最优刚体变换
    (旋转+平移)，返回一个函数，可以把"移动前坐标系"里任意一点映射到"移动后坐标系"。

    用途：PyRosetta对接时只移动蛋白链本身的原子，配体(比如E3配体的exit vector原子)
    如果没有被建成Rosetta的Pose的一部分，不会跟着自动移动——这个函数让我们能用
    "对接前后，蛋白链CA原子的坐标变化"反推出配体应该跟着挪到哪，不需要先把配体
    参数化成Rosetta residue type（那条路线更复杂，见模块docstring的记录）。

    before/after: (N, 3)的numpy数组，N个点在变换前/后的坐标，一一对应。
    """
    before_center, after_center = before.mean(axis=0), after.mean(axis=0)
    before0, after0 = before - before_center, after - after_center
    h = before0.T @ after0
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    correction = np.diag([1.0, 1.0, d])
    rotation = vt.T @ correction @ u.T

    def transform(point: np.ndarray) -> np.ndarray:
        return rotation @ (point - before_center) + after_center

    return transform


def measure_linker_span(linker_smiles: str, n_conformers: int = 50, random_seed: int = 42) -> dict:
    """
    生成n_conformers个构象，测linker两个连接点(裸*标记)之间的真实距离分布——
    用来判断"对接找到的exit vector距离，这根linker够不够长接上"，不是凭经验猜的数字。
    纯RDKit，不需要PyRosetta。
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(linker_smiles)
    dummy_idx = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
    if len(dummy_idx) != 2:
        return {"ok": False, "reason": f"linker应该恰好有2个连接点(*)，实际有{len(dummy_idx)}个"}

    mol_h = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = random_seed
    conf_ids = AllChem.EmbedMultipleConfs(mol_h, numConfs=n_conformers, params=params)
    if len(conf_ids) == 0:
        return {"ok": False, "reason": "构象生成失败"}
    AllChem.MMFFOptimizeMoleculeConfs(mol_h)

    spans = []
    for cid in conf_ids:
        conf = mol_h.GetConformer(cid)
        p1 = np.array(conf.GetAtomPosition(dummy_idx[0]))
        p2 = np.array(conf.GetAtomPosition(dummy_idx[1]))
        spans.append(float(np.linalg.norm(p1 - p2)))

    return {
        "ok": True,
        "n_conformers": len(spans),
        "min_span_angstrom": round(min(spans), 1),
        "max_span_angstrom": round(max(spans), 1),
        "median_span_angstrom": round(float(np.median(spans)), 1),
    }


# ------------------------------------------------------------
# 需要PyRosetta的部分：lazy import，主.venv里诚实返回ok=False
# ------------------------------------------------------------
def run_protein_protein_docking_trials(
    combined_pdb_path: str,
    fixed_chain: str,
    mobile_chain: str,
    mobile_chain_external_point_original: tuple[float, float, float],
    fixed_chain_external_point: tuple[float, float, float],
    n_trials: int = 100,
    random_seed: int | None = None,
) -> dict:
    """
    真实的蛋白-蛋白刚体对接采样：随机朝向+滑入接触，重复n_trials次，每次独立。
    用PyRosetta自带的RigidBodyRandomizeMover + FaDockingSlideIntoContact，
    不需要PatchDock/完整Rosetta C++套件/HPC调度系统。

    mobile_chain_external_point_original: mobile_chain在**原始(未对接)**坐标系里，
    某个不属于蛋白链本身的外部参考点坐标(比如E3配体的exit vector原子)——每次对接后，
    用kabsch_rigid_transform()算出的变换把这个点映射到对接后的坐标系，不需要先把
    配体参数化成Rosetta residue type。
    fixed_chain_external_point: fixed_chain侧的对应参考点，因为fixed_chain在对接
    过程中不移动，这个点全程不变，直接传入即可。

    真实运行过一次（路线B warhead vs 真实CRBN结构4TZ4，n=500）：exit vector最小距离
    26.8Å，0/500次落进任何一个linker的真实span范围内——这是本函数被设计出来之后
    第一次真实暴露出来的发现，不是假设性的，见 doc/routes/route-b-degrader.md。
    """
    try:
        import pyrosetta
        from pyrosetta.rosetta.protocols.docking import FaDockingSlideIntoContact, setup_foldtree
        from pyrosetta.rosetta.protocols.rigid import Partner, RigidBodyRandomizeMover
        from pyrosetta.rosetta.utility import vector1_int
    except ImportError as exc:
        return {
            "ok": False,
            "reason": f"需要PyRosetta做蛋白-蛋白刚体对接，当前解释器未安装: {exc}",
            "combined_pdb_path": combined_pdb_path,
        }

    pyrosetta.init("-mute all")
    base_pose = pyrosetta.pose_from_pdb(combined_pdb_path)
    movable_jumps = vector1_int()
    movable_jumps.append(1)
    setup_foldtree(base_pose, f"{fixed_chain}_{mobile_chain}", movable_jumps)
    scorefxn = pyrosetta.get_fa_scorefxn()

    mobile_begin, mobile_end = base_pose.chain_begin(2), base_pose.chain_end(2)
    original_ca = np.array(
        [
            [base_pose.residue(i).xyz("CA").x, base_pose.residue(i).xyz("CA").y, base_pose.residue(i).xyz("CA").z]
            for i in range(mobile_begin, mobile_end + 1)
        ]
    )
    mobile_point_original = np.array(mobile_chain_external_point_original)
    fixed_point = np.array(fixed_chain_external_point)

    if random_seed is not None:
        pyrosetta.rosetta.numeric.random.rg().set_seed("mt19937", random_seed)

    trials = []
    for trial_idx in range(n_trials):
        pose = pyrosetta.Pose()
        pose.assign(base_pose)
        RigidBodyRandomizeMover(pose, 1, Partner.partner_downstream).apply(pose)
        FaDockingSlideIntoContact(1).apply(pose)

        moved_ca = np.array(
            [
                [pose.residue(i).xyz("CA").x, pose.residue(i).xyz("CA").y, pose.residue(i).xyz("CA").z]
                for i in range(mobile_begin, mobile_end + 1)
            ]
        )
        transform = kabsch_rigid_transform(original_ca, moved_ca)
        mobile_point_new = transform(mobile_point_original)
        exit_vector_distance = float(np.linalg.norm(mobile_point_new - fixed_point))
        trials.append(
            {"trial": trial_idx, "score": round(scorefxn(pose), 1), "exit_vector_distance_angstrom": round(exit_vector_distance, 1)}
        )

    trials.sort(key=lambda t: t["exit_vector_distance_angstrom"])
    return {"ok": True, "n_trials": n_trials, "trials_sorted_by_distance": trials}


def evaluate_ternary_complex_geometry(docking_trials_result: dict, linker_max_span_angstrom: float) -> dict:
    """
    用run_protein_protein_docking_trials()的结果，判断有没有任何一次对接的exit vector
    距离落在某个linker真实能达到的span以内(见measure_linker_span())。
    """
    if not docking_trials_result.get("ok"):
        return docking_trials_result

    trials = docking_trials_result["trials_sorted_by_distance"]
    best = trials[0]
    compatible_trials = [t for t in trials if t["exit_vector_distance_angstrom"] <= linker_max_span_angstrom]
    return {
        "ok": True,
        "n_trials": docking_trials_result["n_trials"],
        "best_exit_vector_distance_angstrom": best["exit_vector_distance_angstrom"],
        "linker_max_span_angstrom": linker_max_span_angstrom,
        "n_geometrically_compatible_trials": len(compatible_trials),
        "any_compatible_pose_found": len(compatible_trials) > 0,
        "note": (
            "找到了几何兼容的姿态，可以把这几个trial的pose拿去做进一步的界面能量精修"
            if compatible_trials
            else f"{docking_trials_result['n_trials']}次独立随机对接里，最好的exit vector距离也有"
            f"{best['exit_vector_distance_angstrom']}Å，超过linker最大span({linker_max_span_angstrom}Å)——"
            "均匀随机刚体采样找不到兼容姿态，不代表三元复合物不存在，更可能是采样策略本身的限制"
            "(真实PRosettaC用PatchDock系统性穷举表面补丁，不是均匀随机)"
        ),
    }


if __name__ == "__main__":
    print("=== run_ternary_complex_stub (预期 ok=False) ===")
    print(run_ternary_complex_stub("demo_warhead", "lenalidomide_n_linked", "peg2"))

    print("\n=== run_protein_protein_docking_trials (主.venv下预期 ok=False，需要.venv310) ===")
    print(
        run_protein_protein_docking_trials(
            combined_pdb_path="routes/route_b_degrader/structures/egfr_crbn_combined.pdb",
            fixed_chain="A",
            mobile_chain="B",
            mobile_chain_external_point_original=(-45.717, 58.320, -98.517),
            fixed_chain_external_point=(-47.554, -0.549, -20.829),
            n_trials=5,
        )
    )

    print("\n=== measure_linker_span (纯RDKit，不需要PyRosetta) ===")
    print("peg2:", measure_linker_span("*CCOCCOCC*"))

    print("\n=== kabsch_rigid_transform 逻辑验证(合成数据：已知的旋转+平移) ===")
    before = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    translation = np.array([5.0, 3.0, 1.0])
    after = before + translation  # 纯平移，没有旋转，最简单的可验证场景
    transform = kabsch_rigid_transform(before, after)
    test_point = np.array([2.0, 2.0, 2.0])
    print(f"测试点{test_point}平移后应该是{test_point + translation}，实际算出来是{transform(test_point)}")

    print("\n=== evaluate_hook_effect_risk 逻辑验证(合成数据) ===")
    print("安全剂量范围(10-100nM，限制性Kd=50nM):")
    print(evaluate_hook_effect_risk(warhead_kd_nm=50, e3_ligand_kd_nm=20, dose_range_nm=(10, 100)))

    print("\n进入风险区间的剂量范围(100-1000nM，限制性Kd=50nM，阈值500nM):")
    print(evaluate_hook_effect_risk(warhead_kd_nm=50, e3_ligand_kd_nm=20, dose_range_nm=(100, 1000)))

    print("\n=== 构造一个完整的 TernaryComplexResult(大部分字段诚实留空) ===")
    result = TernaryComplexResult(
        warhead_id="demo_warhead_from_route_c_scaffold",
        e3_ligand_id="lenalidomide_n_linked",
        linker_id="peg2",
        warhead_kd_nm=50.0,
        e3_ligand_kd_nm=20.0,
    )
    print(result)
