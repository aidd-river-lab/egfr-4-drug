"""
Retrospective验证：用已知答案(真实活性化合物 vs 诱饵/无活性分子)检验打分方法本身
有没有区分能力——这是doc/11-industrial-gap-and-roadmap.md第5节"第一关"要求的那一步，
在本仓库任何对接/打分结果被当真用于排序决策之前，必须先确认这把尺子在公开已知数据上
量得准不准。

两个纯数学指标，不依赖Vina/PyRosetta，可以单独测试：
  - compute_roc_auc(): ROC曲线下面积，衡量"打分能不能把活性分子整体排在诱饵前面"
  - compute_enrichment_factor(): 富集因子，衡量"排名最前面那一小部分里，活性分子
    有没有比随机抽样更集中"——这是虚拟筛选文献里更贴近"实际会怎么用"的指标
    (真实场景里没人会把全部几万个分子都送去合成，只会挑排名前面的一小部分)。

真实数据来源：scripts/run_dude_validation.py跑DUD-E的EGFR数据集产出的
validation/dude_egfr/results/docking_scores.csv，用这两个函数算出最终报告。

第三个函数`compute_auc_by_similarity_to_reference()`是路线A的SHP2验证
(AUC=0.575，明显比EGFR那次弱)跑完之后新增的诊断工具：AUC偏弱时，一个
常见的疑点是"真实数据库里混进了结合在不同位点/不同机制的化合物"(比如
SHP2既有变构tunnel位点抑制剂，也有活性位点抑制剂，强行用同一个受体对接
会制造虚假的低区分度)。这个函数用2D化学结构相似度(Morgan指纹+Tanimoto)
把化合物按"像不像某个已知结合在目标位点的原型分子"分组，分别算AUC——
如果分组后"更像"那组AUC明显更高，说明确实是"混进了别的机制"；如果没有
这个模式(本仓库SHP2那次就是没有)，说明问题更可能是任务本身难，不是
口袋选错了。这是一个比"逐个化合物查文献确认结合位点"便宜得多的初筛，
不能100%取代人工核对，但能先排除/支持这个假设。
"""
from __future__ import annotations

from dataclasses import dataclass


def compute_roc_auc(labels: list[bool], scores: list[float], lower_score_is_more_positive: bool = True) -> dict:
    """
    ROC-AUC，用Mann-Whitney U统计量的等价公式算(不需要sklearn)：
        AUC = U / (n_pos * n_neg)
    其中U是"每一对(正例,负例)里，正例排得比负例更靠前"的次数(并列算半次)。
    这个公式和"曲线下面积"在数学上完全等价，但不需要真的去扫描一条ROC曲线。

    lower_score_is_more_positive=True：对接分数的惯例(Vina/Rosetta都是越负=越强)，
    分数越低越应该被判定为"活性"。如果用的是"分数越高越好"的打分(比如相似度打分)，
    传False。

    labels/scores长度必须相同且一一对应；至少各需要1个正例和1个负例才能算，
    否则diamond退化(AUC没有意义，诚实返回ok=False)。
    """
    if len(labels) != len(scores):
        return {"ok": False, "reason": "labels和scores长度不一致"}
    if len(labels) == 0:
        return {"ok": False, "reason": "没有数据"}

    pos_scores = [s for lbl, s in zip(labels, scores) if lbl]
    neg_scores = [s for lbl, s in zip(labels, scores) if not lbl]
    n_pos, n_neg = len(pos_scores), len(neg_scores)
    if n_pos == 0 or n_neg == 0:
        return {"ok": False, "reason": f"正例({n_pos})和负例({n_neg})必须都至少有1个才能算AUC"}

    # U统计量：对每一对(正例, 负例)，判断正例是否排得比负例更"像活性"(分数更好)。
    # 并列(分数完全相等)算0.5——标准处理，避免系统性偏向任何一边。
    u = 0.0
    for ps in pos_scores:
        for ns in neg_scores:
            if lower_score_is_more_positive:
                better = ps < ns
                tie = ps == ns
            else:
                better = ps > ns
                tie = ps == ns
            if tie:
                u += 0.5
            elif better:
                u += 1.0

    auc = u / (n_pos * n_neg)
    return {
        "ok": True,
        "auc": round(auc, 4),
        "n_pos": n_pos,
        "n_neg": n_neg,
        "interpretation": _interpret_auc(auc),
    }


def _interpret_auc(auc: float) -> str:
    if auc >= 0.9:
        return "优秀区分能力(接近完美排序)"
    if auc >= 0.7:
        return "有意义的区分能力(典型对接方法在公开benchmark上的正常区间)"
    if auc >= 0.55:
        return "弱区分能力，比随机好但不多，排序结果慎用"
    if auc > 0.45:
        return "和随机猜测无法区分，这把尺子在当前数据上没有验证出排序能力"
    return "系统性排反了(AUC显著低于0.5)，检查受体口袋/box坐标/打分方向是不是搞反了"


@dataclass
class EnrichmentResult:
    top_fraction: float
    n_top: int
    n_actives_in_top: int
    n_actives_total: int
    enrichment_factor: float
    max_possible_ef: float  # 如果top里全是活性分子，EF的理论上限(受活性分子总数限制)


def compute_enrichment_factor(
    labels: list[bool], scores: list[float], top_fraction: float, lower_score_is_more_positive: bool = True
) -> dict:
    """
    富集因子(Enrichment Factor)：排名最前面top_fraction比例的分子里，活性分子的
    比例，相对于"随机排序预期比例"的倍数。

        EF = (排名前top_fraction里的活性数 / 排名前top_fraction的总数)
             / (活性总数 / 全部分子总数)

    EF=1表示排序和随机打乱没区别；EF=5表示"如果你只有预算合成测试排名前面的那一批，
    里面活性分子的密度是随机抽样的5倍"——这是虚拟筛选文献里比AUC更贴近真实决策场景
    的指标，因为真实决策从来不是"看完整条排序"，是"只往下走一小段"。
    """
    if len(labels) != len(scores):
        return {"ok": False, "reason": "labels和scores长度不一致"}
    n_total = len(labels)
    if n_total == 0:
        return {"ok": False, "reason": "没有数据"}
    if not (0 < top_fraction <= 1):
        return {"ok": False, "reason": "top_fraction必须在(0, 1]区间"}

    n_actives_total = sum(labels)
    if n_actives_total == 0:
        return {"ok": False, "reason": "没有活性分子(正例)，无法算富集因子"}

    order = sorted(range(n_total), key=lambda i: scores[i], reverse=not lower_score_is_more_positive)
    n_top = max(1, round(n_total * top_fraction))
    top_indices = order[:n_top]
    n_actives_in_top = sum(labels[i] for i in top_indices)

    expected_fraction_actives = n_actives_total / n_total
    observed_fraction_actives = n_actives_in_top / n_top
    ef = observed_fraction_actives / expected_fraction_actives if expected_fraction_actives > 0 else float("nan")
    max_possible_ef = min(n_actives_total, n_top) / n_top / expected_fraction_actives

    result = EnrichmentResult(
        top_fraction=top_fraction,
        n_top=n_top,
        n_actives_in_top=n_actives_in_top,
        n_actives_total=n_actives_total,
        enrichment_factor=round(ef, 2),
        max_possible_ef=round(max_possible_ef, 2),
    )
    return {"ok": True, **result.__dict__}


def compute_auc_by_similarity_to_reference(
    labels: list[bool],
    scores: list[float],
    smiles_list: list[str],
    reference_smiles: str,
    thresholds: list[float] = (0.15, 0.2, 0.25, 0.3),
    lower_score_is_more_positive: bool = True,
) -> dict:
    """
    按"和某个已知参照分子的2D结构相似度"分组，分别算AUC——诊断"AUC整体偏弱，
    是不是因为数据库里混进了结合在别的位点/别的机制的化合物"这个假设。

    reference_smiles：已知确实结合在目标位点的原型分子(比如SHP2的SHP099，
    一定要是查证过真实存在、真实结构的分子，不能是凭记忆编的)。
    thresholds：多个Tanimoto相似度阈值，每个阈值把compounds分成"更像参照分子"
    和"更不像"两组，各自算一次AUC，方便看这个分组是不是真的有区分意义
    (如果每个阈值下两组AUC都差不多，说明相似度分组这个假设没有被数据支持)。

    依赖RDKit(这个仓库其他函数不需要，这里需要算分子指纹)，lazy import。
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, DataStructs
    except ImportError as exc:
        return {"ok": False, "reason": f"需要RDKit算分子相似度: {exc}"}

    ref_mol = Chem.MolFromSmiles(reference_smiles)
    if ref_mol is None:
        return {"ok": False, "reason": f"参照分子SMILES解析失败: {reference_smiles}"}
    ref_fp = AllChem.GetMorganFingerprintAsBitVect(ref_mol, radius=2, nBits=2048)

    if not (len(labels) == len(scores) == len(smiles_list)):
        return {"ok": False, "reason": "labels/scores/smiles_list长度必须一致"}

    similarities = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            similarities.append(None)
            continue
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
        similarities.append(DataStructs.TanimotoSimilarity(ref_fp, fp))

    valid = [(lbl, sc, sim) for lbl, sc, sim in zip(labels, scores, similarities) if sim is not None]
    if not valid:
        return {"ok": False, "reason": "没有能成功算出相似度的化合物"}

    overall_auc = compute_roc_auc(
        [v[0] for v in valid], [v[1] for v in valid], lower_score_is_more_positive
    )

    by_threshold = {}
    for threshold in thresholds:
        high = [(lbl, sc) for lbl, sc, sim in valid if sim >= threshold]
        low = [(lbl, sc) for lbl, sc, sim in valid if sim < threshold]
        by_threshold[threshold] = {
            "high_similarity_n": len(high),
            "high_similarity_auc": compute_roc_auc(
                [h[0] for h in high], [h[1] for h in high], lower_score_is_more_positive
            ),
            "low_similarity_n": len(low),
            "low_similarity_auc": compute_roc_auc(
                [l[0] for l in low], [l[1] for l in low], lower_score_is_more_positive
            ),
        }

    return {"ok": True, "overall_auc": overall_auc, "by_threshold": by_threshold}


if __name__ == "__main__":
    print("=== 完美排序(活性分子分数全部比诱饵低)：AUC应该是1.0 ===")
    labels = [True, True, True, False, False, False]
    scores = [-9.0, -8.5, -8.0, -5.0, -4.5, -4.0]
    print(compute_roc_auc(labels, scores))

    print("\n=== 随机打分(标签和分数完全独立)：n足够大时AUC应该收敛到0.5附近 ===")
    import random

    rng = random.Random(0)
    labels2 = [rng.random() < 0.5 for _ in range(400)]
    scores2 = [rng.uniform(-10, -4) for _ in range(400)]  # 和labels2完全独立生成，没有真实关联
    print(compute_roc_auc(labels2, scores2))

    print("\n=== 排反了(活性分子分数反而更高/更差)：AUC应该接近0.0 ===")
    labels3 = [True, True, True, False, False, False]
    scores3 = [-4.0, -4.5, -5.0, -8.0, -8.5, -9.0]
    print(compute_roc_auc(labels3, scores3))

    print("\n=== 富集因子：前50%刚好全是活性分子 ===")
    print(compute_enrichment_factor(labels, scores, top_fraction=0.5))

    print("\n=== 按相似度分组诊断(真实SHP099 vs 几个虚构的苯环分子做demo) ===")
    shp099 = "CC1(N)CCN(c2cnc(-c3cccc(Cl)c3Cl)c(N)n2)CC1"
    demo_smiles = ["CC1(N)CCN(c2cnc(-c3cccc(Cl)c3Cl)c(N)n2)CC1", "c1ccccc1", "c1ccccc1O", "c1ccccc1N"]
    demo_labels = [True, False, False, True]
    demo_scores = [-9.5, -5.0, -5.2, -8.8]
    print(compute_auc_by_similarity_to_reference(demo_labels, demo_scores, demo_smiles, shp099, thresholds=[0.5]))
