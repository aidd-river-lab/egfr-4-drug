"""
routes/route_a_shp2_sos1/config/rgroup_libraries/*.yaml 的RDKit可解析性验证。
和 test_rgroup_libraries.py(路线C)结构完全对称，独立成文件避免两条路线的
测试互相干扰(呼应 test_rgroup_libraries.py 顶部注释)。
"""
from pathlib import Path

import yaml
from rdkit import Chem

from core.route_config import route_config_dir

CONFIG_DIR = route_config_dir("route_a_shp2_sos1") / "rgroup_libraries"


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
        assert n_dummy == 1, f"{region_name}/{frag_name} 应该恰好有1个连接点(*)，实际有{n_dummy}个: {smiles}"


def test_expected_pocket_regions_present():
    regions = _load_pocket_regions()
    for expected in ["tunnel_amine", "distal_aryl"]:
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
        assert isotopes == [1, 2], f"scaffold {scaffold_name} 应该有同位素标记为1/2的两个连接点，实际: {isotopes}"
