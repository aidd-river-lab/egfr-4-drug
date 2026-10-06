"""
环节4.3真实验证：路线C真实候选(demo_aminopyrimidine_biphenyl骨架)实际上是
reversible_orthosteric(chemistry_route.primary，covalent_warhead_enabled=false)，
不含共价弹头——`core/covalent.py::evaluate_attack_geometry()`从未在真实分子上
跑过，因为真实候选根本没有弹头可测。真正有弹头的只有奥希替尼这个外部参照分子，
而本仓库已有的奥希替尼对接(6LUD/C797S突变体系)用的受体Cys797已经被突变成Ser——
这正是C797S耐药的生物学机制本身，意味着那些受体里压根没有真实的Cys-SG可供测量。

本脚本用真实PDB **4ZAU**(标题原文"AZD9291 COMPLEX WITH WILD TYPE EGFR"，
SEQADV记录只有N端表达标签，790/797位都没有突变记录，Cys797完整保留)：把奥希替尼
当作完整、未反应的小分子非共价对接进这个有完整Cys797的受体，提取对接姿态里
弹头反应碳到Cys797-SG的真实距离和角度，喂给`evaluate_attack_geometry()`——
这是用一个真实世界已知"确实会反应"的正例去检验几何判据函数本身，和DUD-E
retrospective验证是同一个思路(拿已知答案检验方法)。

重要声明：Vina不能模拟共价反应，只给反应前的encounter complex姿态；这里测的是
"弹头有没有摆在合理的攻击轨迹上"，不是真实共价键长(真实共价键长约1.8Å，这里测
的是非共价接触距离)。攻击角的具体原子定义见`core/covalent.py::
extract_michael_acceptor_geometry()`的docstring。

用法：
    .venv310/bin/python scripts/run_covalent_geometry_check.py
"""
from __future__ import annotations

import subprocess
from pathlib import Path

STRUCT_DIR = Path(__file__).resolve().parent.parent / "routes" / "route_c_4th_gen_tki" / "structures"
RAW_PDB = STRUCT_DIR / "4zau_raw.pdb"


def _ligand_centroid_from_pdb(pdb_path: Path, hetcode: str) -> tuple[float, float, float]:
    xs, ys, zs = [], [], []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("HETATM") and line[17:20].strip() == hetcode:
                xs.append(float(line[30:38]))
                ys.append(float(line[38:46]))
                zs.append(float(line[46:54]))
    return (sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs))


def main() -> None:
    from core.structures import clean_chain

    # 1. 清洗：keep_hetatm=True只为了拿YY3(奥希替尼真实共晶坐标)质心做对接盒子中心
    hetatm_path = STRUCT_DIR / "4zau_chainA_with_het.pdb"
    clean_chain(RAW_PDB, chain="A", out_path=hetatm_path, keep_hetatm=True)
    center = _ligand_centroid_from_pdb(hetatm_path, "YY3")
    print(f"[box center] 真实YY3(奥希替尼)共晶坐标质心: {center}")

    # 2. 纯受体(不含HETATM)用于对接，保留完整Cys797
    receptor_pdb = STRUCT_DIR / "4zau_receptor_only.pdb"
    clean_chain(RAW_PDB, chain="A", out_path=receptor_pdb, keep_hetatm=False)

    # 3. meeko准备受体PDBQT
    receptor_pdbqt = STRUCT_DIR / "4zau_receptor.pdbqt"
    result = subprocess.run(
        [
            "mk_prepare_receptor.py",
            "--read_pdb", str(receptor_pdb),
            "-o", str(receptor_pdbqt.with_suffix("")),
            "--allow_bad_res",
            "-p",
        ],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        print("受体准备失败:", result.stderr)
        return
    print(f"[受体准备] 完成: {receptor_pdbqt}")

    # 4. 奥希替尼配体：复用已经验证过的真实SMILES(tests/test_standardize.py验证过InChIKey)
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    from rdkit import Chem
    from rdkit.Chem import AllChem

    osimertinib_smiles = "COc1cc(N(C)CCN(C)C)c(NC(=O)C=C)cc1Nc1nccc(-c2cn(C)c3ccccc23)n1"
    mol = Chem.MolFromSmiles(osimertinib_smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    prep = MoleculePreparation()
    setups = prep.prepare(mol)
    ligand_pdbqt_string, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        print("配体准备失败:", err)
        return
    ligand_pdbqt_path = STRUCT_DIR / "osimertinib_for_4zau.pdbqt"
    ligand_pdbqt_path.write_text(ligand_pdbqt_string)

    # 5. 真实Vina对接
    from vina import Vina

    v = Vina(sf_name="vina", seed=42, verbosity=0)
    v.set_receptor(str(receptor_pdbqt))
    v.compute_vina_maps(center=list(center), box_size=[22.0, 22.0, 22.0])
    v.set_ligand_from_file(str(ligand_pdbqt_path))
    v.dock(exhaustiveness=8, n_poses=1)
    best_affinity = float(v.energies()[0][0])
    docked_path = STRUCT_DIR / "osimertinib_docked_4zau.pdbqt"
    v.write_poses(str(docked_path), n_poses=1, overwrite=True)
    print(f"[对接完成] 最佳打分: {best_affinity} kcal/mol，姿态写入 {docked_path}")

    # 6. 提取真实弹头-Cys797几何
    from core.covalent import CovalentDockingResult, evaluate_attack_geometry, extract_michael_acceptor_geometry

    geometry = extract_michael_acceptor_geometry(
        docked_ligand_pdbqt=docked_path,
        receptor_pdbqt=receptor_pdbqt,
        nucleophile_resname="CYS",
        nucleophile_resnum=797,
        nucleophile_atom_name="SG",
    )
    print("\n[真实几何提取结果]")
    print(geometry)

    if not geometry["ok"]:
        return

    covalent_result = CovalentDockingResult(
        ligand_id="osimertinib_vs_4ZAU_WT_EGFR",
        nucleophile_residue="Cys797",
        warhead_carbon_to_nucleophile_distance_angstrom=geometry["beta_carbon_to_nucleophile_distance_angstrom"],
        attack_angle_degree=geometry["attack_angle_degree"],
        # ki_nm/kinact_per_second/gsh_half_life_hours: 真实值需要湿实验/FEP，本仓库没有，
        # 诚实留空，不编造
    )
    verdict = evaluate_attack_geometry(covalent_result)
    print("\n[evaluate_attack_geometry真实判定]")
    print(verdict)


if __name__ == "__main__":
    main()
