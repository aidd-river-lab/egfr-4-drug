"""
环节6：ADMET 计算描述符 + 结构警示引擎

!! 精度边界声明（请在用这个模块做真实决策前读完）!!
真正的 cLogD7.4 需要"按 pH7.4 下的电离态分布"算分配系数，这需要一个可靠的 pKa 预测器
（ChemAxon/Marvin、或训练好的 pKa 模型）。本模块没有接入 pKa 预测，所以：
  - `clogp` 是 RDKit Crippen.MolLogP，这是中性分子的辛醇/水分配系数，不是 cLogD7.4
  - `clogd_approx` 直接复用 clogp 的数值作为近似占位，会系统性高估带正电荷胺类分子
    在生理pH下的真实分配系数（因为质子化后极性更强、logD应该更低）
  - `pka_proxy` 不是数值pKa预测，只是"是否含有典型强碱性脂肪叔胺"的子结构标记，
    用来粗略估计hERG/CNS渗透的方向性风险，不能替代真实pKa计算
  - `cns_mpo` 因此也是近似值：MW/TPSA/HBD三项是精确计算，cLogP/cLogD两项用同一个
    MolLogP代入、pKa项用上面的粗代理——这会让分数系统性偏乐观或偏保守，具体方向因
    分子而异。真实项目决策前必须换成有pKa支持的工具重新算一遍。
这些近似边界在每个函数的 docstring 里也重复了一遍，是故意的——宁可啰嗦也不要让人
在不知情的情况下把这里的数字当成最终结论。
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import yaml
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

CONFIG_PATH = Path(__file__).resolve().parent.parent / "config" / "structural_alerts.yaml"

# hERG 组合规则(子结构 + 计算属性)用到的"强碱性脂肪叔胺"子结构——不包括芳香胺/酰胺氮/
# 已经被砜化弱化的哌嗪(环节3.2里推荐的那几种弱碱基团不会命中这条)
_BASIC_ALIPHATIC_AMINE_SMARTS = Chem.MolFromSmarts(
    "[NX3;H0;!$(N-[#6]=[O,S,N]);!$(N-a);$(N(-[#6;!a])(-[#6;!a])-[#6;!a])]"
)


def _as_mol(mol_or_smiles) -> Chem.Mol | None:
    if isinstance(mol_or_smiles, Chem.Mol):
        return mol_or_smiles
    return Chem.MolFromSmiles(mol_or_smiles)


# ------------------------------------------------------------------
# 描述符计算
# ------------------------------------------------------------------
def compute_descriptors(mol_or_smiles) -> dict:
    """
    计算环节6.1表格里能用RDKit精确算出来的那几项：MW、clogp(真实是MolLogP不是cLogD)、
    TPSA、HBD、HBA、Fsp3。LLE/LELP 需要外部提供 pic50，不在这里算。
    """
    mol = _as_mol(mol_or_smiles)
    if mol is None:
        return {"ok": False, "error": "无法解析分子"}

    mw = Descriptors.MolWt(mol)
    clogp = Crippen.MolLogP(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    fsp3 = rdMolDescriptors.CalcFractionCSP3(mol)
    has_basic_amine = mol.HasSubstructMatch(_BASIC_ALIPHATIC_AMINE_SMARTS)

    return {
        "ok": True,
        "mw": round(mw, 2),
        "clogp": round(clogp, 2),
        "clogd_approx": round(clogp, 2),  # 近似占位，见模块顶部声明
        "tpsa": round(tpsa, 2),
        "hbd": hbd,
        "hba": hba,
        "fsp3": round(fsp3, 3),
        "has_basic_aliphatic_amine": has_basic_amine,
        "pka_proxy": "high_risk(~8-10)" if has_basic_amine else "low_risk(<7)",
    }


def compute_lle_lelp(pic50: float, descriptors: dict) -> dict:
    """
    LLE = pIC50 - cLogP；LELP = cLogP / LE，LE = 1.37*pIC50 / 重原子数。
    这两个指标必须有真实活性数据(pic50)才有意义，没有活性数据时不要瞎editable造一个。
    """
    if not descriptors.get("ok"):
        return {"ok": False}
    clogp = descriptors["clogp"]
    lle = pic50 - clogp
    return {"ok": True, "lle": round(lle, 2)}


def estimate_herg_risk_proxy(mol_or_smiles) -> dict:
    """
    hERG风险的"子结构+计算属性"组合代理(见 structural_alerts.yaml 底部说明)：
    强碱性脂肪叔胺 + cLogP偏高，是hERG最经典的组合特征。
    这不是hERG IC50预测模型，只是一个方向性的早期预警标记。
    """
    mol = _as_mol(mol_or_smiles)
    if mol is None:
        return {"ok": False}
    desc = compute_descriptors(mol)
    has_amine = desc["has_basic_aliphatic_amine"]
    clogp = desc["clogp"]
    if has_amine and clogp > 3.5:
        risk = "high"
    elif has_amine and clogp > 2.0:
        risk = "medium"
    else:
        risk = "low"
    return {"ok": True, "herg_risk_proxy": risk, "basic_amine": has_amine, "clogp": clogp}


# ------------------------------------------------------------------
# CNS MPO（Pfizer 六参数复合分，Wager et al. 2010/2016）
# 近似边界：见模块顶部声明，cLogP/cLogD用同一个值代入，pKa用粗代理
# ------------------------------------------------------------------
def _desirability_decreasing(x: float, full_score_below: float, zero_score_above: float) -> float:
    if x <= full_score_below:
        return 1.0
    if x >= zero_score_above:
        return 0.0
    return 1.0 - (x - full_score_below) / (zero_score_above - full_score_below)


def _desirability_range(x: float, low_full: float, low_edge: float, high_full: float, high_edge: float) -> float:
    if low_full <= x <= high_full:
        return 1.0
    if x < low_full:
        return 0.0 if x <= low_edge else (x - low_edge) / (low_full - low_edge)
    return 0.0 if x >= high_edge else 1.0 - (x - high_full) / (high_edge - high_full)


def compute_cns_mpo(descriptors: dict) -> dict:
    if not descriptors.get("ok"):
        return {"ok": False}
    mw = descriptors["mw"]
    clogp = descriptors["clogp"]
    clogd = descriptors["clogd_approx"]
    tpsa = descriptors["tpsa"]
    hbd = descriptors["hbd"]
    pka_high_risk = descriptors["has_basic_aliphatic_amine"]

    s_mw = _desirability_decreasing(mw, full_score_below=360, zero_score_above=500)
    s_clogp = _desirability_decreasing(clogp, full_score_below=3, zero_score_above=5)
    s_clogd = _desirability_decreasing(clogd, full_score_below=2, zero_score_above=4)
    s_tpsa = _desirability_range(tpsa, low_full=40, low_edge=20, high_full=90, high_edge=120)
    s_hbd = _desirability_decreasing(hbd, full_score_below=1, zero_score_above=3)
    s_pka = 0.5 if pka_high_risk else 1.0  # 粗代理：有强碱性叔胺打5折，不是真实pKa desirability曲线

    total = s_mw + s_clogp + s_clogd + s_tpsa + s_hbd + s_pka
    return {
        "ok": True,
        "cns_mpo": round(total, 2),
        "components": {
            "mw": round(s_mw, 2),
            "clogp": round(s_clogp, 2),
            "clogd_approx": round(s_clogd, 2),
            "tpsa": round(s_tpsa, 2),
            "hbd": round(s_hbd, 2),
            "pka_proxy": round(s_pka, 2),
        },
        "pass_threshold_4_0": total >= 4.0,
        "caveat": "clogd/pka为近似值，见core/admet.py模块顶部声明",
    }


# ------------------------------------------------------------------
# 结构警示引擎
# ------------------------------------------------------------------
@dataclass
class AlertHit:
    rule_name: str
    severity: str
    action: str
    advice: str


@dataclass
class AlertReport:
    hits: list[AlertHit] = field(default_factory=list)

    @property
    def verdict(self) -> str:
        """三态：任意一条REJECT就整体REJECT；否则任意一条OPTIMIZE_NEEDED就整体OPTIMIZE_NEEDED；都没有就PASS。"""
        actions = {h.action for h in self.hits}
        if "REJECT" in actions:
            return "REJECT"
        if "OPTIMIZE_NEEDED" in actions:
            return "OPTIMIZE_NEEDED"
        return "PASS"


def _load_alert_rules() -> dict:
    with open(CONFIG_PATH, encoding="utf-8") as f:
        return yaml.safe_load(f)["rules"]


_ALERT_RULES_CACHE: dict | None = None


def run_structural_alerts(mol_or_smiles, warhead_smarts: str | None = None) -> AlertReport:
    """
    跑 structural_alerts.yaml 里的全部规则，返回三态判定(PASS/OPTIMIZE_NEEDED/REJECT)。

    warhead_smarts: 如果候选分子走共价路线(route B/C)，设计好的弹头本身会命中
    michael_acceptor_offtarget 规则——这是预期的，不应该被当成缺陷。传入弹头的SMARTS，
    命中这个SMARTS的迈克尔受体会被豁免，不进最终报告。不要为了消掉这个警示去改规则本身
    （见 structural_alerts.yaml 里 michael_acceptor_offtarget 的 advice 字段）。
    """
    global _ALERT_RULES_CACHE
    mol = _as_mol(mol_or_smiles)
    if mol is None:
        return AlertReport(hits=[AlertHit("parse_error", "high", "REJECT", "分子无法解析")])

    if _ALERT_RULES_CACHE is None:
        _ALERT_RULES_CACHE = _load_alert_rules()
    rules = _ALERT_RULES_CACHE

    warhead_patt = Chem.MolFromSmarts(warhead_smarts) if warhead_smarts else None

    hits = []
    for name, rule in rules.items():
        base_patt = Chem.MolFromSmarts(rule["smarts"])
        if base_patt is None or not mol.HasSubstructMatch(base_patt):
            continue
        # 共价弹头豁免：仅对迈克尔受体规则生效
        if name == "michael_acceptor_offtarget" and warhead_patt is not None:
            if mol.HasSubstructMatch(warhead_patt):
                continue
        exempted = False
        for ex_smarts in rule.get("exempt_if", []) or []:
            ex_patt = Chem.MolFromSmarts(ex_smarts)
            if ex_patt is not None and mol.HasSubstructMatch(ex_patt):
                exempted = True
                break
        if exempted:
            continue
        hits.append(AlertHit(rule_name=name, severity=rule["severity"], action=rule["action"], advice=rule["advice"]))

    return AlertReport(hits=hits)


if __name__ == "__main__":
    OSIMERTINIB = "COc1cc(N(C)CCN(C)C)c(NC(=O)C=C)cc1Nc1nccc(-c2cn(C)c3ccccc23)n1"
    mol = Chem.MolFromSmiles(OSIMERTINIB)

    desc = compute_descriptors(mol)
    print("=== descriptors (osimertinib) ===")
    for k, v in desc.items():
        print(f"  {k}: {v}")

    print("\n=== CNS MPO ===")
    mpo = compute_cns_mpo(desc)
    print(mpo)

    print("\n=== hERG risk proxy ===")
    print(estimate_herg_risk_proxy(mol))

    print("\n=== structural alerts (不传warhead，弹头会被标记) ===")
    report = run_structural_alerts(mol)
    for h in report.hits:
        print(f"  [{h.severity}] {h.rule_name} -> {h.action}: {h.advice}")
    print(f"  verdict: {report.verdict}")

    print("\n=== structural alerts (传入丙烯酰胺弹头SMARTS做豁免) ===")
    report2 = run_structural_alerts(mol, warhead_smarts="[CX3]=[CX3][CX3]=[OX1]")
    for h in report2.hits:
        print(f"  [{h.severity}] {h.rule_name} -> {h.action}: {h.advice}")
    print(f"  verdict: {report2.verdict}")
