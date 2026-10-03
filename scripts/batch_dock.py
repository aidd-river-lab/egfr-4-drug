"""
批量跑真实Vina对接：给一组compound_id+SMILES、一个receptor pdbqt、一个口袋中心，
生成3D构象->PDBQT->对接->落库funnel_scores。必须用 .venv310 跑(需要真实meeko/vina)。

用法：
    .venv310/bin/python scripts/batch_dock.py <candidates.json> <receptor.pdbqt> <ensemble_id> <structures_dir> <cx> <cy> <cz>
"""
import json
import sys

sys.path.insert(0, ".")

from rdkit import Chem
from rdkit.Chem import AllChem
import subprocess
from core.docking import run_vina_docking
from db.connect import get_connection
from db.repository import upsert_funnel_score


def main():
    candidates_path, receptor_pdbqt, ensemble_id, out_dir, cx, cy, cz = sys.argv[1:8]
    with open(candidates_path) as f:
        candidates = json.load(f)
    center = (float(cx), float(cy), float(cz))

    conn = get_connection()
    for cid, smi in candidates.items():
        mol = Chem.MolFromSmiles(smi)
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        if AllChem.EmbedMolecule(mol, params) != 0:
            print(f"{cid}: 3D embed失败，跳过")
            continue
        AllChem.MMFFOptimizeMolecule(mol)
        sdf_path = f"{out_dir}/{cid}_3d.sdf"
        writer = Chem.SDWriter(sdf_path)
        writer.write(mol)
        writer.close()

        pdbqt_path = f"{out_dir}/{cid}.pdbqt"
        proc = subprocess.run(
            [".venv310/bin/mk_prepare_ligand.py", "-i", sdf_path, "-o", pdbqt_path],
            capture_output=True, text=True,
        )
        if proc.returncode != 0:
            print(f"{cid}: ligand prep失败: {proc.stderr}")
            continue

        out_pdbqt = f"{out_dir}/{cid}_docked.pdbqt"
        result = run_vina_docking(
            receptor_pdbqt=receptor_pdbqt, ligand_pdbqt=pdbqt_path, center=center,
            box_size=(20.0, 20.0, 20.0), out_pdbqt=out_pdbqt, exhaustiveness=8, n_poses=3,
        )
        print(f"{cid} -> {result.best_affinity_kcal_mol} kcal/mol, ok={result.ok}")
        if result.ok:
            upsert_funnel_score(
                conn, compound_id=cid, ensemble_id=ensemble_id, funnel_level="L1",
                method="vina_docking", score_value=result.best_affinity_kcal_mol,
                error_estimate_kcal_mol=2.5, n_poses_or_frames=result.n_poses,
                passed_filter=result.best_affinity_kcal_mol <= -6.0, raw_output_path=out_pdbqt,
            )
    conn.close()


if __name__ == "__main__":
    main()
