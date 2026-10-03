"""
config/rgroup_libraries/*.yaml 里每个SMILES片段的RDKit可解析性验证。
这是 pocket_regions.yaml 顶部注释承诺的那份测试脚本。

重点检查两件事(呼应 core/enumerate.py 里molzip那个曾经静默失败的bug)：
  1. 每个片段都能被RDKit解析+sanitize
  2. 每个片段恰好有一个虚拟原子(*)作为连接点——多于或少于一个都会让molzip的结果不可控
"""
from pathlib import Path

import yaml
from rdkit import Chem

CONFIG_DIR = Path(__file__).resolve().parent.parent / "config" / "rgroup_libraries"


def _load_pocket_regions():
    with open(CONFIG_DIR / "pocket_regions.yaml", encoding="utf-8") as f:
        return yaml.safe_load(f)


def _load_scaffolds():
    with open(CONFIG_DIR / "scaffolds.yaml", encoding="utf-8") as f:
        return yaml.safe_load(f)


def _all_fragments():
    regions = _load_pocket_regions()
    for region_name, region in regions.items():
        if not isinstance(region, dict) or "fragments" not in region:
            continue
        for frag in region["fragments"]:
            yield region_name, frag["name"], frag["smiles"]


def test_all_fragments_parse():
    for region_name, frag_name, smiles in _all_fragments():
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None, f"{region_name}/{frag_name} 的SMILES无法解析: {smiles}"


def test_all_fragments_have_exactly_one_dummy_atom():
    for region_name, frag_name, smiles in _all_fragments():
        mol = Chem.MolFromSmiles(smiles)
        n_dummy = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 0)
        assert n_dummy == 1, (
            f"{region_name}/{frag_name} 应该恰好有1个连接点(*)，实际有{n_dummy}个: {smiles}"
        )


def test_at_least_one_region_per_expected_pocket():
    regions = _load_pocket_regions()
    for expected in ["hinge", "entrance_797", "hydrophobic_back_pocket", "solvent_exposed"]:
        assert expected in regions
        assert len(regions[expected]["fragments"]) >= 1


def test_scaffold_core_smiles_parses_with_isotope_labeled_dummies():
    scaffolds = _load_scaffolds()
    for scaffold_name, scaffold in scaffolds.items():
        if not isinstance(scaffold, dict) or "core_smiles" not in scaffold:
            continue
        mol = Chem.MolFromSmiles(scaffold["core_smiles"])
        assert mol is not None, f"scaffold {scaffold_name} 的core_smiles无法解析"
        isotopes = sorted(atom.GetIsotope() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 0)
        assert isotopes == [1, 2, 3], (
            f"scaffold {scaffold_name} 应该有同位素标记为1/2/3的三个连接点，实际: {isotopes}"
        )
