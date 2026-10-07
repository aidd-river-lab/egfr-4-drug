"""
路线A的第二块retrospective验证：SOS1(之前只做了SHP2那一半，这次补齐)。

数据来源：ChEMBL真实查询，target_chembl_id=CHEMBL2079846(人源SOS1，已核实
organism=Homo sapiens)。拉取标准化IC50/Ki(standard_relation='='，单位nM)
的全部活性记录，按molecule_chembl_id去重(中位数)，过滤data_validity_comment
标记有问题的记录(1条/1252条)。

真实类别分布，老实记录一个和SHP2不一样的发现：active(IC50<=1μM)有776个，
但inactive(>=10μM)只有**19个**——负例远比SHP2那次(279个)少。这不是抽样
问题，是真实数据本身的分布：SOS1抑制剂是比SHP2更新、更聚焦(主要是2018年
后KRAS相关研究驱动)的领域，历史上"做出来但失败"的化合物远没有SHP2那么多
被详细记录。这意味着这次算出来的AUC，负例这一侧的统计功效天然更弱，解读
结果时要把这个考虑进去，不能因为n_neg小就提高SHP2的可信度对比标准。

active随机抽样(种子42)到400个，inactive全部保留(19个)，总共419个真实
化合物，对接进真实6SCM受体(BI-3406结合的催化邻近口袋，box中心取自
真实晶体坐标质心，这次会话已经验证过能重新对接出BI-3406真实姿态)。

用法：
    .venv310/bin/python scripts/run_chembl_sos1_validation.py
"""
from __future__ import annotations

import csv
import json
import time
from pathlib import Path

LABELED_SET_PATH = Path(__file__).resolve().parent.parent / "validation" / "chembl_sos1" / "raw" / "labeled_set.json"
RESULTS_DIR = Path(__file__).resolve().parent.parent / "validation" / "chembl_sos1" / "results"
RECEPTOR_PDBQT = (
    Path(__file__).resolve().parent.parent / "routes" / "route_a_shp2_sos1" / "structures" / "6scm_receptor.pdbqt"
)
# 真实BI-3406(配体L7H)在6SCM里的质心坐标，这次会话已经验证过能重新对接出真实姿态
BOX_CENTER = (1.01, -31.62, -44.43)
BOX_SIZE = (20.0, 20.0, 20.0)
RANDOM_SEED = 42
EXHAUSTIVENESS = 8


def _prepare_ligand_pdbqt(smiles: str) -> str | None:
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=RANDOM_SEED, useRandomCoords=True) != 0:
        return None
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception:
        pass
    try:
        setups = MoleculePreparation().prepare(mol)
        pdbqt_string, ok, _err = PDBQTWriterLegacy.write_string(setups[0])
    except Exception:
        return None
    return pdbqt_string if ok else None


def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    with open(LABELED_SET_PATH) as f:
        labeled = json.load(f)

    compounds = [(mid, v["smiles"], v["label"], v["median_nm"]) for mid, v in labeled.items()]
    print(f"[数据] {sum(1 for _,_,l,_ in compounds if l=='active')}个真实active + "
          f"{sum(1 for _,_,l,_ in compounds if l=='inactive')}个真实inactive(均来自ChEMBL真实测量)")

    from vina import Vina

    v = Vina(sf_name="vina", seed=RANDOM_SEED, verbosity=0)
    v.set_receptor(str(RECEPTOR_PDBQT))
    t0 = time.time()
    v.compute_vina_maps(center=list(BOX_CENTER), box_size=list(BOX_SIZE))
    print(f"[受体网格] 计算完成，{time.time() - t0:.1f}秒")

    out_path = RESULTS_DIR / "docking_scores.csv"
    n_ok, n_fail = 0, 0
    t_start = time.time()
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["molecule_chembl_id", "label", "median_ic50_or_ki_nm", "smiles", "best_affinity_kcal_mol", "ok", "error"])
        for i, (mid, smiles, label, median_nm) in enumerate(compounds):
            pdbqt_string = None
            error = ""
            try:
                pdbqt_string = _prepare_ligand_pdbqt(smiles)
                if pdbqt_string is None:
                    error = "ligand_prep_failed"
            except Exception as exc:
                error = f"ligand_prep_exception: {exc}"

            best_affinity = None
            if pdbqt_string is not None:
                try:
                    v.set_ligand_from_string(pdbqt_string)
                    v.dock(exhaustiveness=EXHAUSTIVENESS, n_poses=1)
                    best_affinity = float(v.energies()[0][0])
                except Exception as exc:
                    error = f"docking_exception: {exc}"

            ok = best_affinity is not None
            n_ok += int(ok)
            n_fail += int(not ok)
            writer.writerow([mid, label, median_nm, smiles, best_affinity, ok, error])

            if (i + 1) % 50 == 0:
                f.flush()
                elapsed = time.time() - t_start
                rate = (i + 1) / elapsed
                remaining = (len(compounds) - i - 1) / rate if rate > 0 else float("nan")
                print(f"[进度] {i + 1}/{len(compounds)} (成功{n_ok}/失败{n_fail}) "
                      f"已用{elapsed / 60:.1f}分钟，预计剩余{remaining / 60:.1f}分钟")

    print(f"\n[完成] 成功{n_ok} / 失败{n_fail}，结果写入 {out_path}")


if __name__ == "__main__":
    main()
