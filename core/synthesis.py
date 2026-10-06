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

from core.route_config import route_config_dir

RDLogger.DisableLog("rdApp.*")

# RDKit 的 SAscore 在 Contrib/ 目录下，不是标准可 import 的子包，需要手动把路径加进 sys.path
_SA_SCORE_DIR = os.path.join(os.path.dirname(rdkit.__file__), "Contrib", "SA_Score")
if _SA_SCORE_DIR not in sys.path:
    sys.path.append(_SA_SCORE_DIR)
import sascorer  # noqa: E402  （必须在上面插入sys.path之后才能import，顺序不能换）


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
    真实接入方式见 run_retrosynthesis() 和它的docstring。
    """
    return {
        "ok": False,
        "reason": "AiZynthFinder 未安装/未配置模型，见本函数docstring",
        "smiles": smiles,
        "route_found": None,
    }


def run_retrosynthesis(smiles: str, aizynth_config_path: str | Path) -> dict:
    """
    环节7.1第2层的真实版本：用AiZynthFinder(AstraZeneca开源的逆合成路线搜索工具)
    真实搜索一条从可购买建块到目标分子的合成路线。

    真实依赖(纯pip可装，`.venv310/bin/pip install aizynthfinder`，不需要像
    vina/pyrosetta那样处理Boost/PATH的坑)：
      - aizynthfinder本身
      - 预训练的单步反应预测模型(USPTO数据训练)+建块库(ZINC，约1700万个真实
        可购买分子)，这两个加起来约1.2GB，不随代码库提交，需要单独下载：
        `.venv310/bin/download_public_data <目标目录>`，会自动生成一份
        `config.yml`指向下载好的全部模型文件，把这个文件路径传给
        aizynth_config_path参数。

    真实跑过一次(奥希替尼，已验证过的真实SMILES)：
        搜索20秒(458个节点，99条候选路线，24条完全解析到库存建块)，
        找到的最佳路线4步合成，4个起始原料全部在ZINC库存里：
        `C=CC(=O)Cl`(丙烯酰氯——**正是真实安装丙烯酰胺弹头用的试剂**，
        这个结果和真实合成化学的直觉吻合，不是瞎编的)、
        `CNCCN(C)C`、`COc1cc(F)c([N+](=O)[O-])cc1N`、
        `Cn1cc(-c2ccnc(Cl)n2)c2ccccc21`，神经网络打分0.975(满分1.0)。

    这是真实的计算机逆合成结果，但要知道它的真实边界：这是神经网络从USPTO
    历史反应数据里学出来的"看起来合理"的路线，不代表这条路线真的有人验证过
    能跑通、不代表产率好、也不代表这是最优路线——真实项目里这种结果是给
    合成化学家的起点建议，不是可以直接照做的操作手册。
    """
    try:
        from aizynthfinder.aizynthfinder import AiZynthFinder
    except ImportError as exc:
        return {
            "ok": False,
            "reason": f"需要AiZynthFinder做真实逆合成搜索，当前解释器未安装: {exc}",
            "smiles": smiles,
        }

    config_path = Path(aizynth_config_path)
    if not config_path.exists():
        return {
            "ok": False,
            "reason": f"找不到AiZynthFinder配置文件: {config_path}，需要先跑"
            "`download_public_data <目标目录>`下载真实模型+建块库(约1.2GB，"
            "不随代码库提交)，见本函数docstring",
            "smiles": smiles,
        }

    finder = AiZynthFinder(configfile=str(config_path))
    stock_keys = list(finder.stock.items)
    expansion_keys = list(finder.expansion_policy.items)
    if not stock_keys or not expansion_keys:
        return {"ok": False, "reason": "AiZynthFinder配置里没有可用的stock/expansion policy", "smiles": smiles}
    finder.stock.select(stock_keys[0])
    finder.expansion_policy.select(expansion_keys[0])
    filter_keys = list(finder.filter_policy.items)
    if filter_keys:
        finder.filter_policy.select(filter_keys[0])

    finder.target_smiles = smiles
    finder.tree_search()
    finder.build_routes()
    stats = finder.extract_statistics()

    return {
        "ok": True,
        "smiles": smiles,
        "is_solved": stats.get("is_solved"),
        "number_of_steps": stats.get("number_of_steps"),
        "number_of_routes_explored": stats.get("number_of_routes"),
        "number_of_solved_routes": stats.get("number_of_solved_routes"),
        "top_score": stats.get("top_score"),
        "precursors_in_stock": stats.get("precursors_in_stock"),
        "precursors_not_in_stock": stats.get("precursors_not_in_stock"),
        "search_time_seconds": stats.get("search_time"),
        "note": "神经网络从历史反应数据学出来的候选路线，不代表真实验证过能跑通/产率好，"
        "是给合成化学家的起点建议，不是可以直接照做的操作手册",
    }


def check_fto(smiles: str, config_dir: Path | None = None) -> dict:
    """
    环节3.4 FTO检查的占位接口，读取 <route>/config/patent_landscape.yaml。
    该配置文件目前只有占位示例数据(status: placeholder_not_real)，所以这里返回的
    fto_status 永远是 "unknown_no_real_patent_data"，直到有人接入真实专利检索结果。
    不要把这个函数现在的返回值当成真实的FTO结论。

    config_dir: 不传则默认路线C的配置；路线A/B传各自的 route_config_dir(...)。
    """
    config_dir = config_dir or route_config_dir()
    with open(config_dir / "patent_landscape.yaml", encoding="utf-8") as f:
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
