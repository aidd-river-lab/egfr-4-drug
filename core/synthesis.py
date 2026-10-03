"""
环节7：合成可及性

三层筛查里，本模块真正实现的是第1层（启发式分数，秒级，RDKit contrib SAscore）。
第2层（AiZynthFinder计算机逆合成）和FTO的Markush匹配需要额外的模型/数据，本模块
提供正确的函数签名和"未就绪"的明确报错，不用随机数或假路线冒充结果——
参见 config/patent_landscape.yaml 顶部的声明，这是同一个原则的延续。
"""
from __future__ import annotations

import os
import sys
from dataclasses import dataclass
from pathlib import Path

import rdkit
import yaml
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

# RDKit 的 SAscore 在 Contrib/ 目录下，不是标准可 import 的子包，需要手动把路径加进 sys.path
_SA_SCORE_DIR = os.path.join(os.path.dirname(rdkit.__file__), "Contrib", "SA_Score")
if _SA_SCORE_DIR not in sys.path:
    sys.path.append(_SA_SCORE_DIR)
import sascorer  # noqa: E402  （必须在上面插入sys.path之后才能import，顺序不能换）

CONFIG_DIR = Path(__file__).resolve().parent.parent / "config"


def compute_sa_score(mol_or_smiles) -> dict:
    """
    RDKit SAscore：1(容易合成)~10(很难合成)的启发式分数，秒级，只能用来排除明显离谱的，
    不能当决策依据(环节7.1)。
    """
    mol = mol_or_smiles if isinstance(mol_or_smiles, Chem.Mol) else Chem.MolFromSmiles(mol_or_smiles)
    if mol is None:
        return {"ok": False, "error": "无法解析分子"}
    score = sascorer.calculateScore(mol)
    return {"ok": True, "sa_score": round(score, 2), "interpretation": "1=容易, 10=很难"}


@dataclass
class SynthesisRecord:
    """环节7.3 要求必须记录的字段。这些字段大多需要逆合成工具或药化人工填写，
    本dataclass只定义schema，不是每个字段都有自动计算方法——没有数据就留None，
    不要为了"看起来完整"去瞎填默认值。
    """

    compound_id: str
    smiles: str
    sa_score: float | None = None
    n_steps: int | None = None  # 需要AiZynthFinder或药化人工路线规划
    longest_linear_sequence: int | None = None
    key_bb_cas: str | None = None
    bb_lead_time_days: int | None = None
    estimated_cost_per_20mg: float | None = None
    chiral_resolution_needed: bool | None = None
    scale_up_risk: str | None = None  # 自由文本，比如"含叠氮步骤"
    fto_status: str | None = None  # in_claim / near_claim / clear，见 patent_landscape.yaml


def run_retrosynthesis_stub(smiles: str) -> dict:
    """
    环节7.1第2层：计算机逆合成(AiZynthFinder)的占位接口。

    本仓库的开发环境没有安装 AiZynthFinder（需要额外下载预训练的单步反应模型，
    体积较大，且需要和具体的可购建块库配套），所以这里不冒充一个假路线返回给调用方。
    真实接入方式：
        pip install aizynthfinder
        下载官方预训练模型 + 建块库(见 AiZynthFinder 文档)
        用 aizynthfinder.aizynthfinder.AiZynthFinder 类替换这个函数体
    """
    return {
        "ok": False,
        "reason": "AiZynthFinder 未安装/未配置模型，见本函数docstring",
        "smiles": smiles,
        "route_found": None,
    }


def check_fto(smiles: str) -> dict:
    """
    环节3.4 FTO检查的占位接口，读取 config/patent_landscape.yaml。
    该配置文件目前只有占位示例数据(status: placeholder_not_real)，所以这里返回的
    fto_status 永远是 "unknown_no_real_patent_data"，直到有人接入真实专利检索结果。
    不要把这个函数现在的返回值当成真实的FTO结论。
    """
    with open(CONFIG_DIR / "patent_landscape.yaml", encoding="utf-8") as f:
        landscape = yaml.safe_load(f)

    real_patents = [p for p in landscape["competitor_patents"] if p.get("status") != "placeholder_not_real"]
    if not real_patents:
        return {
            "ok": True,
            "fto_status": "unknown_no_real_patent_data",
            "smiles": smiles,
            "note": "patent_landscape.yaml 里还没有真实专利数据，见该文件顶部声明",
        }
    # 真实数据接入后，这里才应该做 Markush SMARTS 匹配
    return {"ok": True, "fto_status": "not_implemented_pending_real_data", "smiles": smiles}


if __name__ == "__main__":
    osimertinib = "COc1cc(N(C)CCN(C)C)c(NC(=O)C=C)cc1Nc1nccc(-c2cn(C)c3ccccc23)n1"
    print("=== SA score ===")
    print(compute_sa_score(osimertinib))

    print("\n=== retrosynthesis stub (预期: ok=False, 诚实报告未就绪) ===")
    print(run_retrosynthesis_stub(osimertinib))

    print("\n=== FTO check (预期: unknown_no_real_patent_data) ===")
    print(check_fto(osimertinib))

    print("\n=== SynthesisRecord schema demo ===")
    rec = SynthesisRecord(
        compound_id="DEMO-001",
        smiles=osimertinib,
        sa_score=compute_sa_score(osimertinib)["sa_score"],
    )
    print(rec)
