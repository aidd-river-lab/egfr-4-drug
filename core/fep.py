"""
环节4 L4：FEP 相对结合自由能

真实计算需要 OpenFE/OpenMM 或 Schrödinger FEP+，单对微扰 2-8 GPU·小时，
50对的map需要200-400 GPU·小时(见doc里的算力预算表)，本环境无GPU。

这个模块里真实、可测试的部分是 FEPMapResult.qc_pass —— "用实验锚点验证整张map
是否可信"的质控逻辑(设计文档环节4.4："预测值和实测值的MUE > 1.5 kcal/mol时，
整张map作废，回去查位姿或质子化态")，这条规则不需要真的跑过FEP就能写对、测对。
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import yaml
from pathlib import Path

PIPELINE_CONFIG = Path(__file__).resolve().parent.parent / "config" / "pipeline.yaml"


@dataclass
class FEPPerturbation:
    ligand_a: str
    ligand_b: str
    predicted_ddg_kcal_mol: float
    uncertainty_kcal_mol: float
    experimental_ddg_kcal_mol: float | None = None  # 非None表示这对是实验锚点


@dataclass
class FEPMapResult:
    genotype: str
    perturbations: list[FEPPerturbation] = field(default_factory=list)

    @property
    def anchor_perturbations(self) -> list[FEPPerturbation]:
        return [p for p in self.perturbations if p.experimental_ddg_kcal_mol is not None]

    def mue_vs_experiment(self) -> float | None:
        anchors = self.anchor_perturbations
        if not anchors:
            return None
        errors = [abs(p.predicted_ddg_kcal_mol - p.experimental_ddg_kcal_mol) for p in anchors]
        return float(np.mean(errors))

    def qc_pass(self) -> dict:
        """
        设计文档环节4.4的硬规则：没有实验锚点的map不能用；MUE超过阈值整张map作废。
        这是真实的质控逻辑，用合成数据就能完整验证正确性（见 __main__）。
        """
        with open(PIPELINE_CONFIG, encoding="utf-8") as f:
            threshold = yaml.safe_load(f)["funnel"]["L4"]["qc"]["mue_reject_threshold_kcal_mol"]

        if not self.anchor_perturbations:
            return {"ok": False, "reason": "没有实验锚点，这张FEP map不能用来做决策(环节4.4硬要求)"}
        mue = self.mue_vs_experiment()
        if mue > threshold:
            return {"ok": False, "reason": f"MUE={mue:.2f} kcal/mol 超过阈值{threshold}，整张map作废，回去查位姿或质子化态", "mue": mue}
        return {"ok": True, "mue": mue, "n_anchors": len(self.anchor_perturbations)}


def run_fep_stub(genotype: str, perturbation_pairs: list[tuple[str, str]]) -> dict:
    """占位：真实运行需要OpenFE/OpenMM或FEP+，本环境无GPU未执行，不编造ΔΔG数值。"""
    return {
        "ok": False,
        "reason": "需要OpenFE/OpenMM(或商业FEP+) + GPU集群，本环境无GPU未执行",
        "genotype": genotype,
        "n_pairs_requested": len(perturbation_pairs),
    }


if __name__ == "__main__":
    print("=== run_fep_stub (预期 ok=False) ===")
    print(run_fep_stub("del19_C797S", [("A", "B"), ("B", "C")]))

    print("\n=== qc_pass 逻辑验证(合成数据) ===")
    good_map = FEPMapResult(
        genotype="del19_C797S",
        perturbations=[
            FEPPerturbation("A", "B", predicted_ddg_kcal_mol=-0.8, uncertainty_kcal_mol=0.3, experimental_ddg_kcal_mol=-1.0),
            FEPPerturbation("B", "C", predicted_ddg_kcal_mol=0.4, uncertainty_kcal_mol=0.2, experimental_ddg_kcal_mol=0.6),
            FEPPerturbation("C", "D", predicted_ddg_kcal_mol=-1.2, uncertainty_kcal_mol=0.4),  # 非锚点
        ],
    )
    print("good_map:", good_map.qc_pass())

    bad_map = FEPMapResult(
        genotype="del19_C797S",
        perturbations=[
            FEPPerturbation("A", "B", predicted_ddg_kcal_mol=2.5, uncertainty_kcal_mol=0.3, experimental_ddg_kcal_mol=-1.0),
        ],
    )
    print("bad_map:", bad_map.qc_pass())

    no_anchor_map = FEPMapResult(
        genotype="del19_C797S",
        perturbations=[FEPPerturbation("A", "B", predicted_ddg_kcal_mol=-0.5, uncertainty_kcal_mol=0.3)],
    )
    print("no_anchor_map:", no_anchor_map.qc_pass())
