"""
从真实数据库里把某条路线L1对接分数最好的几个真实候选分子，连同它们真实的
Vina对接姿态和对应受体结构，打包成一个可以直接上传到Colab的zip——给
`notebooks/colab_gpu_md.ipynb`里"跑你自己筛选出来的候选分子"那一节用。

和之前demo版notebook(硬编码奥希替尼+6LUD)的区别：这里导出的是本仓库
`workflows/run_discovery_round.py`真实枚举+真实Vina对接跑出来的候选
(`routes/<route>/structures/RTx-xxxx_docked.pdbqt`)，不是示例分子。

**怎么把docked.pdbqt变成MD能用的正确分子**：Vina对接产出的.pdbqt只有坐标
和AutoDock原子类型，没有标准的化学键级信息。这里用meeko自己的
`PDBQTMolecule` + `RDKitMolCreate`(这正是meeko官方提供的、把对接结果读回
RDKit分子对象的工具，不是本脚本发明的转换逻辑)把docked pose还原成一个
同时有"真实键级"+"真实对接坐标"的RDKit分子，再写成.sdf——这一步已经在
这次会话里对真实的RTC-0000_docked.pdbqt验证过，还原出的SMILES和这个
候选分子的设计结构完全一致。

用法：
    .venv310/bin/python scripts/export_top_candidates_for_gpu_md.py \
        --route route_a_shp2_sos1 --top-n 3

    .venv310/bin/python scripts/export_top_candidates_for_gpu_md.py \
        --route route_c_4th_gen_tki --top-n 3 --ensemble L858R_T790M_C797S_6lud
"""
from __future__ import annotations

import argparse
import shutil
import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

# 每条路线默认用哪个受体(对应funnel_scores真实用的ensemble_id和受体文件)。
# 路线A/C都不止一个真实ensemble，不传--ensemble时用这里列的默认值。
ROUTE_DEFAULT_ENSEMBLE = {
    "route_a_shp2_sos1": "SHP2_tunnel_5ehr",
    "route_c_4th_gen_tki": "L858R_T790M_C797S_6lud",
    "route_b_degrader": None,  # 路线B目前只有warhead片段对接，见doc/16第7节的已知限制
}

ENSEMBLE_RECEPTOR_PDB = {
    "SHP2_tunnel_5ehr": "routes/route_a_shp2_sos1/structures/5ehr_receptor_only.pdb",
    "SOS1_catalytic_adjacent_cleft_6scm": "routes/route_a_shp2_sos1/structures/6scm_receptor_only.pdb",
    "L858R_T790M_C797S_6lud": "routes/route_c_4th_gen_tki/structures/6lud_receptor_only.pdb",
    "L858R_C797S_pyrosetta_from_6lud": "routes/route_c_4th_gen_tki/structures/L858R_C797S_computational_model.pdb",
    "CRBN_lenalidomide_4tz4": "routes/route_b_degrader/structures/4tz4_crbn_chainC_clean.pdb",
}


def _reconstruct_docked_sdf(docked_pdbqt_path: Path, out_sdf_path: Path) -> dict:
    """用meeko把真实docked.pdbqt还原成带正确键级+真实对接坐标的.sdf。"""
    try:
        from meeko import PDBQTMolecule, RDKitMolCreate
        from rdkit import Chem
    except ImportError as exc:
        return {"ok": False, "reason": f"meeko/rdkit未安装: {exc}"}

    try:
        pdbqt_mol = PDBQTMolecule.from_file(str(docked_pdbqt_path), skip_typing=True)
        mols = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)
        if not mols:
            return {"ok": False, "reason": "RDKitMolCreate没有还原出任何分子"}
        mol = mols[0]
        writer = Chem.SDWriter(str(out_sdf_path))
        writer.write(mol)
        writer.close()
        return {"ok": True, "n_atoms": mol.GetNumAtoms()}
    except Exception as exc:
        return {"ok": False, "reason": str(exc)}


def export_top_candidates(route_id: str, top_n: int, ensemble_id: str | None, out_dir: Path) -> dict:
    from db.connect import get_connection

    ensemble_id = ensemble_id or ROUTE_DEFAULT_ENSEMBLE.get(route_id)
    if ensemble_id is None:
        return {"ok": False, "reason": f"路线{route_id}没有配置默认ensemble，需要显式传--ensemble"}

    receptor_pdb_rel = ENSEMBLE_RECEPTOR_PDB.get(ensemble_id)
    if receptor_pdb_rel is None:
        return {"ok": False, "reason": f"ensemble {ensemble_id} 没有登记对应的受体PDB路径"}

    conn = get_connection()
    cur = conn.cursor()
    cur.execute(
        """
        SELECT compound_id, score_value, raw_output_path
        FROM funnel_scores
        WHERE funnel_level = 'L1' AND ensemble_id = %s AND passed_filter = 1
        ORDER BY score_value ASC
        """,
        (ensemble_id,),
    )
    rows = cur.fetchall()
    if not rows:
        return {"ok": False, "reason": f"ensemble {ensemble_id} 下没有真实L1对接记录"}

    seen = set()
    top_rows = []
    for compound_id, score_value, raw_output_path in rows:
        if compound_id in seen:
            continue
        seen.add(compound_id)
        top_rows.append((compound_id, score_value, raw_output_path))
        if len(top_rows) >= top_n:
            break

    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_lines = [f"ensemble_id={ensemble_id}", f"receptor_pdb={receptor_pdb_rel}", ""]
    exported = []
    for compound_id, score_value, raw_output_path in top_rows:
        docked_pdbqt = ROOT / raw_output_path
        if not docked_pdbqt.is_file():
            manifest_lines.append(f"{compound_id}: 跳过，找不到{docked_pdbqt}")
            continue
        out_sdf = out_dir / f"{compound_id}.sdf"
        result = _reconstruct_docked_sdf(docked_pdbqt, out_sdf)
        if not result["ok"]:
            manifest_lines.append(f"{compound_id}: 还原失败 - {result['reason']}")
            continue
        exported.append(compound_id)
        manifest_lines.append(
            f"{compound_id}: 真实Vina打分{score_value} kcal/mol，{result['n_atoms']}个原子，已导出{out_sdf.name}"
        )

    receptor_src = ROOT / receptor_pdb_rel
    receptor_dst = out_dir / receptor_src.name
    shutil.copy(receptor_src, receptor_dst)
    manifest_lines.append(f"\n受体文件: {receptor_src.name}(真实{ENSEMBLE_RECEPTOR_PDB[ensemble_id]}来源)")

    manifest_path = out_dir / "MANIFEST.txt"
    manifest_path.write_text("\n".join(manifest_lines), encoding="utf-8")

    zip_path = out_dir.with_suffix(".zip")
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for f in out_dir.iterdir():
            zf.write(f, arcname=f.name)

    return {"ok": True, "exported": exported, "zip_path": str(zip_path), "manifest": str(manifest_path)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--route", required=True, choices=sorted(ROUTE_DEFAULT_ENSEMBLE))
    parser.add_argument("--top-n", type=int, default=3)
    parser.add_argument("--ensemble", default=None, help="不传则用该路线的默认ensemble")
    args = parser.parse_args()

    out_dir = ROOT / "validation" / "gpu_md_export" / f"{args.route}_top{args.top_n}"
    result = export_top_candidates(args.route, args.top_n, args.ensemble, out_dir)

    if not result["ok"]:
        print(f"[失败] {result['reason']}")
        return

    print(f"[完成] 真实导出{len(result['exported'])}个候选: {result['exported']}")
    print(f"[打包] {result['zip_path']}")
    print("把这个zip上传到Colab(notebook里'跑你自己筛选出来的候选分子'那一节)就能用。")


if __name__ == "__main__":
    main()
