"""
路线A专属retrospective验证：DUD-E不收录SHP2/SOS1(已核对完整102个靶点列表，
见doc/15)，改用ChEMBL真实活性数据自己构造一个标签数据集，跑同样的
AUC/富集因子检验。

数据来源：ChEMBL REST API，真实查询(非离线缓存)，target_chembl_id=CHEMBL3864
(人源PTPN11/SHP2，已核实organism=Homo sapiens，和5EHR里的SHP2蛋白一致)。
拉取标准化IC50/Ki(standard_relation='='，单位nM)的全部活性记录，按
molecule_chembl_id去重(多条测量取中位数)，过滤掉data_validity_comment标记
为有问题的记录(27条/3453条)。

真实类别分布(这是SHP2/ChEMBL数据和DUD-E很不一样的地方，老实记录)：
按IC50<=1000nM(1μM)算active、>=10000nM(10μM)算inactive，筛出1420个active、
只有279个inactive——**活性分子远多于无活性分子，和DUD-E刻意构造的"多诱饵少
活性"完全相反**，这是药物化学文献本身的发表偏向(大部分被登记的化合物是某个
优化项目里相对成功的，"做出来但不行"的化合物远比"做出来且有效"的少被详细
报道/登记)，不是本仓库的抽样偏差。为了控制对接总耗时，active随机抽样
(种子42)到400个，inactive全部保留(279个，这是更稀缺的一类，没有理由再抽样
减少)，总计679个真实化合物。

**和DUD-E验证的关键方法学差异，必须讲清楚**：DUD-E的"诱饵"是计算挑选出来、
物理性质相似但结构故意不同的分子，用来模拟"随机遇到的、大概率不结合"的
阴性对照；这里的inactive是**真实测过、来自同一批医药化学项目的弱效/无效
化合物**，很可能和active在化学结构上高度相似(同一个优化系列里的类似物)——
这是一个比DUD-E**更难**的区分任务(结构相似的分子更难靠对接分数区分开)，
不是更容易，解读AUC数字时要考虑这一点。

另一个必须声明的局限：5EHR代表的是SHP2的**变构tunnel位点**(SHP099结合的
位点)，但ChEMBL里"SHP2抑制剂"这个大类历史上还包含一类完全不同机制的
**催化活性位点抑制剂**(直接堵催化半胱氨酸，和tunnel位点变构抑制剂是两类
不同的化学空间/结合模式)。本脚本没有逐个化合物去查文献确认绑定位点，如果
labeled_set里混入了催化位点抑制剂，它们对tunnel位点对接打分不好完全可能是
"结合模式本来就不对"，不是"这个方法没有区分力"——这会让AUC看起来比真实
情况更差，解读时要考虑这个混淆因素。

用法：
    .venv310/bin/python scripts/run_chembl_shp2_validation.py
"""
from __future__ import annotations

import csv
import json
import time
from pathlib import Path

LABELED_SET_PATH = Path(__file__).resolve().parent.parent / "validation" / "chembl_shp2" / "raw" / "labeled_set.json"
RESULTS_DIR = Path(__file__).resolve().parent.parent / "validation" / "chembl_shp2" / "results"
RECEPTOR_PDBQT = (
    Path(__file__).resolve().parent.parent / "routes" / "route_a_shp2_sos1" / "structures" / "5ehr_receptor.pdbqt"
)
# 真实SHP099(配体5OD)在5EHR里的质心坐标，从5ehr_chainA_clean.pdb直接量出来的，
# 不是凭空设的；box_size沿用本仓库route A之前对接16个真实候选时用的20x20x20。
BOX_CENTER = (22.42, 41.19, 4.70)
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
