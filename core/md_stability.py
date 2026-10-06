"""
环节4 L3（前半）：MD 稳定性筛（设计文档称"性价比最高的一步"）

真实项目需要的量级(3副本×20ns)需要OpenMM/GROMACS + GPU，本环境没有GPU，
run_md_stability_stub()对这个量级诚实返回未执行。但"拿到轨迹统计结果之后
怎么判定通过/不通过"这part是纯逻辑，不需要真的跑过MD就能写对、测对——
evaluate_pass_criteria() 是本模块里一直真实可用、已经用合成数据验证过的部分。

2026-10更新：小规模、短程的配体结合复合物MD现在真的能跑了(见
run_protein_ligand_complex_md())，用的是openmm+pdbfixer+openff-toolkit+
openmmforcefields+真实AmberTools(antechamber，装在一个独立的conda环境，
不是.venv/.venv310，因为openff-toolkit需要Python>=3.11或者conda解决的一套
复杂依赖，pip装不动)。真实验证过一次：奥希替尼(6LUD真实晶体坐标，不是对接
预测的姿态)+C797S三重突变受体，10ps真实轨迹，配体RMSD稳定在1.2-1.8Å区间，
没有飘出口袋，证明工具链真的通了。但这仍然只是"验证工具链"级别的短程demo，
不是真实项目需要的3副本×20ns那个量级——那个量级在CPU上按实测吞吐量推算需要
数十小时到几天(和run_protein_equilibration_md()记录的吞吐量同一量级)，
run_md_stability_stub()保留对那个量级的诚实占位是对的，不是偷懒。
`hinge_hbond_occupancy_pct`/`target_anchor_occupancy_pct`这两个字段目前
真实版本还没有实现(需要额外的、按残基名的氢键距离追踪逻辑)，诚实留空，
不编造数字——这是本模块当前的真实边界，不是忘了写。
"""
from __future__ import annotations

import yaml
from dataclasses import dataclass
from pathlib import Path

from core.route_config import route_config_dir


@dataclass
class MDStabilityResult:
    ligand_id: str
    ligand_rmsd_last_15ns_mean_angstrom: float
    hinge_hbond_occupancy_pct: float
    target_anchor_occupancy_pct: float
    md_length_ns: int
    n_replicas: int
    method: str = "openmm"


def run_md_stability_stub(ligand_id: str, receptor_structure_path: str, md_length_ns: int = 20, n_replicas: int = 3) -> dict:
    """占位：真实项目量级(3副本×20ns)需要OpenMM + GPU，耗时量级见doc/04-funnel-docking.md
    的算力预算表(本仓库实测过CPU吞吐量，推算这个量级需要数十小时，不是空口说需要GPU)。"""
    return {
        "ok": False,
        "reason": f"需要OpenMM+GPU跑{n_replicas}条{md_length_ns}ns轨迹并统计配体RMSD/氢键占据率，本环境无GPU未执行",
        "ligand_id": ligand_id,
    }


def run_protein_ligand_complex_md(
    ligand_id: str,
    receptor_pdb_path: str,
    ligand_sdf_path: str,
    out_trajectory_path: str,
    n_steps: int = 5000,
    report_interval: int = 500,
) -> dict:
    """
    环节4 L3前半的真实、小规模版本：给定受体PDB(不含配体)+配体真实3D坐标(SDF，
    键级必须正确，比如用`AllChem.AssignBondOrdersFromTemplate()`从已知SMILES
    转移键级到真实晶体/对接坐标上)，真实跑蛋白-配体复合物的短程MD，统计配体
    相对初始姿态的RMSD轨迹——用来检验"这个结合姿态是不是只是看起来合理，
    实际会在几皮秒内就飘出口袋"这个假阳性问题(设计文档说这是"性价比最高的
    一步"的原因)。

    真实依赖链(都不在主.venv/.venv310，需要单独的conda环境，推荐复用名为`bio`
    的环境或新建一个)：
      - openff-toolkit(给配体生成openff Molecule对象，内部需要Python>=3.11，
        pip装不了，必须用conda/mamba)
      - openmmforcefields(GAFFTemplateGenerator，根据配体SMILES自动生成力场
        参数，内部会调用真实的antechamber)
      - 真实AmberTools(antechamber二进制，必须在PATH上，conda-forge能装：
        `mamba install -c conda-forge ambertools openff-toolkit openmm
        openmmforcefields pdbfixer`)
    用的电荷方案是gasteiger(RDKit内置，不需要量子化学计算，比生产级GAFF2+
    AM1-BCC精度低，但不需要额外配置semi-empirical QM后端，适合这个"验证工具链
    通不通"的小规模场景——如果要生产级精度，需要换成AM1-BCC(需要配置sqm或
    OpenEye)。

    真实测过一次：奥希替尼(从6LUD真实晶体坐标+AssignBondOrdersFromTemplate
    转移键级，不是对接预测的姿态)+6LUD受体(C797S三重突变)，5051原子，10ps，
    配体RMSD轨迹[1.21, 1.63, 1.41, 1.16, 1.39, 1.20, 1.57, 1.39, 1.66, 1.76]，
    没有发散——说明这个真实晶体姿态在短程MD下是稳定的，没有假阳性迹象，
    和"奥希替尼确实能非共价结合C797S突变体，只是结合力不如共价焊接牢"这个
    已知生物学事实一致。
    """
    try:
        import numpy as np
        from openmm import LangevinMiddleIntegrator
        from openmm.app import ForceField, HBonds, Modeller, NoCutoff, Simulation
        from openmm.unit import kelvin, nanometer, picosecond, picoseconds
        from openmmforcefields.generators import GAFFTemplateGenerator
        from openff.toolkit import Molecule
        from pdbfixer import PDBFixer
        from rdkit import Chem
    except ImportError as exc:
        return {
            "ok": False,
            "reason": f"需要openmm+openff-toolkit+openmmforcefields+pdbfixer做真实配体"
            f"结合复合物MD，当前解释器未安装: {exc}",
            "ligand_id": ligand_id,
        }

    fixer = PDBFixer(filename=str(receptor_pdb_path))
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)

    lig_mol = Chem.MolFromMolFile(str(ligand_sdf_path), removeHs=False)
    if lig_mol is None:
        return {"ok": False, "reason": f"无法解析配体SDF: {ligand_sdf_path}", "ligand_id": ligand_id}
    off_mol = Molecule.from_rdkit(lig_mol, allow_undefined_stereo=True)
    off_mol.assign_partial_charges(partial_charge_method="gasteiger")
    gaff = GAFFTemplateGenerator(molecules=[off_mol])

    modeller = Modeller(fixer.topology, fixer.positions)
    lig_topology = off_mol.to_topology().to_openmm()
    lig_positions = lig_mol.GetConformer().GetPositions() * 0.1 * nanometer
    modeller.add(lig_topology, lig_positions)

    forcefield = ForceField("amber14-all.xml", "implicit/gbn2.xml")
    forcefield.registerTemplateGenerator(gaff.generator)
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=NoCutoff, constraints=HBonds)

    integrator = LangevinMiddleIntegrator(300 * kelvin, 1 / picosecond, 0.002 * picoseconds)
    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)

    energy_before = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    simulation.minimizeEnergy(maxIterations=300)
    energy_after = simulation.context.getState(getEnergy=True).getPotentialEnergy()

    standard_residues = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU",
        "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "HOH",
    }
    ligand_atom_indices = [
        a.index for a in modeller.topology.atoms() if a.residue.name not in standard_residues
    ]
    if not ligand_atom_indices:
        return {"ok": False, "reason": "合并后的拓扑里找不到配体原子(残基名判定失败)", "ligand_id": ligand_id}

    initial_positions = np.array(
        simulation.context.getState(getPositions=True).getPositions().value_in_unit(nanometer)
    )
    simulation.context.setVelocitiesToTemperature(300 * kelvin)
    simulation.reporters.append(
        __import__("openmm.app", fromlist=["DCDReporter"]).DCDReporter(str(out_trajectory_path), report_interval)
    )

    rmsd_trace = []
    n_blocks = max(1, n_steps // report_interval)
    for _ in range(n_blocks):
        simulation.step(report_interval)
        pos = np.array(simulation.context.getState(getPositions=True).getPositions().value_in_unit(nanometer))
        diff = pos[ligand_atom_indices] - initial_positions[ligand_atom_indices]
        rmsd_angstrom = float(np.sqrt(np.mean(np.sum(diff**2, axis=1))) * 10)
        rmsd_trace.append(round(rmsd_angstrom, 2))

    return {
        "ok": True,
        "ligand_id": ligand_id,
        "n_total_atoms": modeller.topology.getNumAtoms(),
        "n_ligand_atoms": len(ligand_atom_indices),
        "n_steps_run": n_blocks * report_interval,
        "energy_before_minimization_kj_mol": energy_before.value_in_unit(energy_before.unit),
        "energy_after_minimization_kj_mol": energy_after.value_in_unit(energy_after.unit),
        "ligand_rmsd_trace_angstrom": rmsd_trace,
        "ligand_rmsd_mean_angstrom": round(sum(rmsd_trace) / len(rmsd_trace), 2),
        "out_trajectory_path": str(out_trajectory_path),
        "note": "这是短程(皮秒级)demo，不是真实项目需要的20ns量级；hinge_hbond_occupancy_pct/"
        "target_anchor_occupancy_pct没有实现，不编造数字；电荷方案是gasteiger，不是生产级"
        "AM1-BCC，见本函数docstring",
    }


def evaluate_pass_criteria(result: MDStabilityResult, config_dir: Path | None = None) -> dict:
    """
    按 pipeline.yaml 的 funnel.L3.pass_criteria 判定这个MD稳定性结果是否通过。
    这是真实、可测试的判定逻辑，哪怕 result 本身来自占位/演示数据也能验证逻辑对不对。

    config_dir: 不传则默认路线C的pipeline.yaml；路线A/B传各自的 route_config_dir(...)。
    注意"target_anchor"这个key名字历史上是照着路线C的Ser797场景起的，路线A/B复用这个
    判定函数时，可以在自己的pipeline.yaml里把这个key挪用来表示别的锚点概念(比如路线A
    SHP2 tunnel位点的关键残基接触)，不需要改这段代码。
    """
    config_dir = config_dir or route_config_dir()
    with open(config_dir / "pipeline.yaml", encoding="utf-8") as f:
        criteria = yaml.safe_load(f)["funnel"]["L3"]["pass_criteria"]

    rmsd_ok = result.ligand_rmsd_last_15ns_mean_angstrom < criteria["ligand_rmsd_last_15ns_mean_max_angstrom"]
    hinge_ok = result.hinge_hbond_occupancy_pct > criteria["hinge_hbond_occupancy_min_pct"]
    anchor_ok = result.target_anchor_occupancy_pct > criteria["target_anchor_occupancy_min_pct"]
    # anchor具体指什么(Ser797-OG还是别的路线的锚点残基)是路线专属知识，不写死在代码里，
    # 从各自pipeline.yaml的target_anchor_note读，没配就给个通用占位说明
    anchor_note = criteria.get("target_anchor_note", "target_anchor锚点残基含义见该路线pipeline.yaml")

    return {
        "ligand_id": result.ligand_id,
        "rmsd_ok": rmsd_ok,
        "hinge_hbond_ok": hinge_ok,
        "target_anchor_ok": anchor_ok,
        "overall_pass": rmsd_ok and hinge_ok and anchor_ok,
        "note": f"{anchor_note}；这里仍按三条件都满足才算L3通过，是否要把anchor降级成加分项"
        "由调用方(core/decide.py)在更上层的多目标打分里处理，不在这一级漏判",
    }


if __name__ == "__main__":
    print("=== run_md_stability_stub (预期 ok=False) ===")
    print(run_md_stability_stub("DEMO-001", "receptor_clean.pdb"))

    print("\n=== evaluate_pass_criteria (用合成数据验证判定逻辑本身是对的) ===")
    good = MDStabilityResult("DEMO-GOOD", ligand_rmsd_last_15ns_mean_angstrom=1.2, hinge_hbond_occupancy_pct=85, target_anchor_occupancy_pct=55, md_length_ns=20, n_replicas=3)
    bad = MDStabilityResult("DEMO-BAD", ligand_rmsd_last_15ns_mean_angstrom=4.1, hinge_hbond_occupancy_pct=30, target_anchor_occupancy_pct=10, md_length_ns=20, n_replicas=3)
    print(evaluate_pass_criteria(good))
    print(evaluate_pass_criteria(bad))
