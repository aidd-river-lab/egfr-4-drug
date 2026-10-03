"""
路线B核心难点：warhead-EGFR-E3连接酶三元复合物预测

二元对接(warhead结合EGFR、E3配体结合CRBN)只是必要条件，不是充分条件——
真正决定降解效率的是三元复合物能不能形成，以及形成时warhead/E3配体的结合界面
会不会互相干扰(负协同)还是互相促进(正协同)。这是二元对接完全捕捉不到的，
需要PRosettaC/Rosetta这类专门做蛋白-蛋白-小分子三元复合物对接的工具，本环境
未安装，run_ternary_complex_stub()诚实占位。

本模块真实、可测试的部分：hook effect(钩状效应)判定。这是双功能降解剂公认的
真实药理学现象——高浓度下warhead和E3配体各自的二元结合各自趋于饱和，
反而会让"同一个蛋白同时结合两边"的三元复合物比例下降，剂量-响应曲线呈钟形
而不是单调递增。这条判据是纯浓度/平衡常数的数学关系，不需要真的跑过三元复合物
3D预测就能验证对不对（设计文档外的通用PROTAC药理学原理，见
doc/routes/route-b-degrader.md 的引用）。
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass
class TernaryComplexResult:
    """三元复合物评估结果的字段schema。真实3D预测部分留空(None)，不编造数值。"""

    warhead_id: str
    e3_ligand_id: str
    linker_id: str
    warhead_kd_nm: float | None = None  # warhead对EGFR的二元结合常数(可以从L1对接/FEP拿到)
    e3_ligand_kd_nm: float | None = None  # E3配体对CRBN的二元结合常数(通常是已知文献值)
    cooperativity_alpha: float | None = None  # 协同因子：>1正协同，<1负协同，=1无协同；需要真实三元复合物数据才能测
    predicted_dc50_nm: float | None = None
    predicted_dmax_pct: float | None = None


def evaluate_hook_effect_risk(
    warhead_kd_nm: float, e3_ligand_kd_nm: float, dose_range_nm: tuple[float, float]
) -> dict:
    """
    Hook effect判定：当给药浓度上限超过两个二元结合常数中较大者的若干倍时，
    两个二元结合会分别趋于饱和，三元复合物浓度反而下降(钟形剂量-响应曲线)。
    这里用一个保守的经验阈值(10倍)标记进入风险区间的起点，不是精确预测IC50会跌多少——
    精确预测需要完整的三元复合物平衡方程求解(还需要未知的cooperativity_alpha)，
    这里只做方向性的风险标记，诚实反映计算能力的边界。

    warhead_kd_nm/e3_ligand_kd_nm: 两个二元结合的平衡解离常数
    dose_range_nm: 计划考察的给药/实测浓度范围 (low, high)
    """
    limiting_kd = max(warhead_kd_nm, e3_ligand_kd_nm)
    hook_risk_threshold_nm = 10 * limiting_kd
    dose_low, dose_high = dose_range_nm

    in_risk_zone = dose_high >= hook_risk_threshold_nm
    safe_margin_fold = hook_risk_threshold_nm / dose_low if dose_low > 0 else None

    return {
        "limiting_kd_nm": limiting_kd,
        "hook_risk_threshold_nm": hook_risk_threshold_nm,
        "dose_range_nm": dose_range_nm,
        "hook_effect_risk": in_risk_zone,
        "safe_margin_fold": round(safe_margin_fold, 1) if safe_margin_fold else None,
        "note": (
            "剂量上限超过两个二元结合常数中较大者的10倍，进入hook effect风险区间"
            if in_risk_zone
            else "剂量范围在风险阈值以内，但这只是方向性判断，不替代真实细胞DC50/Dmax曲线"
        ),
    }


def run_ternary_complex_stub(warhead_id: str, e3_ligand_id: str, linker_id: str) -> dict:
    """占位：真实三元复合物3D预测需要PRosettaC/Rosetta类工具，本环境未安装，不编造cooperativity数值。"""
    return {
        "ok": False,
        "reason": "需要PRosettaC/Rosetta做蛋白-蛋白-小分子三元复合物对接，本环境未安装，不编造cooperativity_alpha/DC50数值",
        "warhead_id": warhead_id,
        "e3_ligand_id": e3_ligand_id,
        "linker_id": linker_id,
    }


if __name__ == "__main__":
    print("=== run_ternary_complex_stub (预期 ok=False) ===")
    print(run_ternary_complex_stub("demo_warhead", "lenalidomide_n_linked", "peg2"))

    print("\n=== evaluate_hook_effect_risk 逻辑验证(合成数据) ===")
    print("安全剂量范围(10-100nM，限制性Kd=50nM):")
    print(evaluate_hook_effect_risk(warhead_kd_nm=50, e3_ligand_kd_nm=20, dose_range_nm=(10, 100)))

    print("\n进入风险区间的剂量范围(100-1000nM，限制性Kd=50nM，阈值500nM):")
    print(evaluate_hook_effect_risk(warhead_kd_nm=50, e3_ligand_kd_nm=20, dose_range_nm=(100, 1000)))

    print("\n=== 构造一个完整的 TernaryComplexResult(大部分字段诚实留空) ===")
    result = TernaryComplexResult(
        warhead_id="demo_warhead_from_route_c_scaffold",
        e3_ligand_id="lenalidomide_n_linked",
        linker_id="peg2",
        warhead_kd_nm=50.0,
        e3_ligand_kd_nm=20.0,
    )
    print(result)
