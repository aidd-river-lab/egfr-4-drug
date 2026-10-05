"""
Retrospective验证：DUD-E的EGFR数据集(542个真实已知活性化合物 + 35050个诱饵分子，
来自dude.docking.org, egfr target, 受体基于PDB 2RGP)，检验core/docking.py的
Vina+meeko流程本身有没有区分能力——这是doc/11-industrial-gap-and-roadmap.md
第5节建议的"第一关"，必须先跑这个再决定要不要砸更多算力/数据。

数据来源：validation/dude_egfr/raw/，真实从dude.docking.org下载
(2012年发表的DUD-E数据集，引用见Mysinger et al. J. Med. Chem. 2012)。
receptor.pdb是DUD-E官方提供的、不带element列的老式PDB格式，RDKit能宽松解析，
meeko的严格解析器不能——先用RDKit读入再重新导出成标准格式(receptor_fixed.pdb)，
再喂给mk_prepare_receptor.py（已经做过一次，结果在raw/receptor.pdbqt）。

box中心/大小：从crystal_ligand.mol2算出的配体质心(16.79, 34.53, 91.29)，
盒子22x22x22Å(配体本身延展9-14Å，留够余量)。

decoys数量太大(35050个)，单机对接不现实(按测出来的~3.15秒/配体，全跑完要30+
小时)，这里固定随机种子抽样2000个(约3.7:1的诱饵:活性比例，远低于DUD-E原生
50:1，但542个真阳性样本本身已经足够撑起一个有统计功效的AUC估计，抽样决策
诚实记录在这里，不是偷偷减少样本量)。

用法：
    .venv310/bin/python scripts/run_dude_validation.py
跑完后用 core/validation.py 的 compute_roc_auc / compute_enrichment_factor
从输出的 docking_scores.csv 算最终指标。
"""
from __future__ import annotations

import csv
import random
import time
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

RAW_DIR = Path(__file__).resolve().parent.parent / "validation" / "dude_egfr" / "raw"
RESULTS_DIR = Path(__file__).resolve().parent.parent / "validation" / "dude_egfr" / "results"
RECEPTOR_PDBQT = RAW_DIR / "receptor.pdbqt"
BOX_CENTER = (16.79, 34.53, 91.29)  # crystal_ligand.mol2质心，见本文件docstring
BOX_SIZE = (22.0, 22.0, 22.0)
N_DECOYS_SAMPLE = 2000
RANDOM_SEED = 42
EXHAUSTIVENESS = 8  # 和本仓库其余对接调用保持一致(core/docking.py的默认值)


def _parse_ism(path: Path) -> list[tuple[str, str]]:
    """DUD-E的.ism文件格式：每行"SMILES zinc_id/编号 [chembl_id]"，空格分隔。返回[(smiles, id), ...]。"""
    rows = []
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 2:
                continue
            rows.append((parts[0], parts[1]))
    return rows


def _prepare_ligand_pdbqt(smiles: str) -> str | None:
    """从SMILES生成3D构象+Gasteiger电荷+PDBQT字符串。失败返回None，不抛异常中断整批。"""
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    embed_result = AllChem.EmbedMolecule(mol, randomSeed=RANDOM_SEED, useRandomCoords=True)
    if embed_result != 0:
        return None
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception:
        pass  # 优化失败不致命，直接用embed的构象

    prep = MoleculePreparation()
    try:
        setups = prep.prepare(mol)
        pdbqt_string, ok, _err = PDBQTWriterLegacy.write_string(setups[0])
    except Exception:
        return None
    return pdbqt_string if ok else None


def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    actives = _parse_ism(RAW_DIR / "actives_final.ism")
    all_decoys = _parse_ism(RAW_DIR / "decoys_final.ism")
    rng = random.Random(RANDOM_SEED)
    decoys = rng.sample(all_decoys, min(N_DECOYS_SAMPLE, len(all_decoys)))

    print(f"[数据] {len(actives)}个真实活性化合物(全量) + {len(decoys)}个诱饵(从{len(all_decoys)}个里随机抽样)")

    from vina import Vina

    v = Vina(sf_name="vina", seed=RANDOM_SEED, verbosity=0)
    v.set_receptor(str(RECEPTOR_PDBQT))
    t0 = time.time()
    v.compute_vina_maps(center=list(BOX_CENTER), box_size=list(BOX_SIZE))
    print(f"[受体网格] 计算完成，{time.time() - t0:.1f}秒，整个批次复用同一份网格(不逐配体重算)")

    out_path = RESULTS_DIR / "docking_scores.csv"
    compounds = [(smi, cid, "active") for smi, cid in actives] + [(smi, cid, "decoy") for smi, cid in decoys]

    n_ok, n_fail = 0, 0
    t_start = time.time()
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["compound_id", "source", "smiles", "best_affinity_kcal_mol", "ok", "error"])
        for i, (smiles, compound_id, source) in enumerate(compounds):
            pdbqt_string = None
            error = ""
            try:
                pdbqt_string = _prepare_ligand_pdbqt(smiles)
                if pdbqt_string is None:
                    error = "ligand_prep_failed(embed/meeko)"
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
            writer.writerow([compound_id, source, smiles, best_affinity, ok, error])

            if (i + 1) % 50 == 0:
                f.flush()
                elapsed = time.time() - t_start
                rate = (i + 1) / elapsed
                remaining = (len(compounds) - i - 1) / rate if rate > 0 else float("nan")
                print(
                    f"[进度] {i + 1}/{len(compounds)} (成功{n_ok}/失败{n_fail})"
                    f" 已用{elapsed / 60:.1f}分钟，预计剩余{remaining / 60:.1f}分钟"
                )

    print(f"\n[完成] 成功{n_ok} / 失败{n_fail}，结果写入 {out_path}")


if __name__ == "__main__":
    main()
