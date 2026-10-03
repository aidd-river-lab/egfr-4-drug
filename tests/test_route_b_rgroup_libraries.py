"""
routes/route_b_degrader/config/rgroup_libraries/*.yaml 的RDKit可解析性验证。
和路线A/C的pocket_regions验证同一类测试，但路线B有三类片段(各自连接点数量不同)：
warhead=1个连接点，linker=2个连接点，e3_ligand=1个连接点。
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

from core.enumerate import load_bifunctional_fragments
from core.route_config import route_config_dir

CONFIG_DIR = route_config_dir("route_b_degrader")


def test_warheads_parse_with_exactly_one_dummy():
    libs = load_bifunctional_fragments(CONFIG_DIR)
    for frag in libs["warheads"]:
        mol = Chem.MolFromSmiles(frag["smiles"])
        assert mol is not None, f"warhead {frag['name']} 无法解析"
        n_dummy = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 0)
        assert n_dummy == 1, f"warhead {frag['name']} 应该有1个连接点，实际{n_dummy}"


def test_linkers_parse_with_exactly_two_dummies():
    libs = load_bifunctional_fragments(CONFIG_DIR)
    for frag in libs["linkers"]:
        mol = Chem.MolFromSmiles(frag["smiles"])
        assert mol is not None, f"linker {frag['name']} 无法解析"
        n_dummy = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 0)
        assert n_dummy == 2, f"linker {frag['name']} 应该有2个连接点，实际{n_dummy}"


def test_e3_ligands_parse_with_exactly_one_dummy():
    libs = load_bifunctional_fragments(CONFIG_DIR)
    for frag in libs["e3_ligands"]:
        mol = Chem.MolFromSmiles(frag["smiles"])
        assert mol is not None, f"e3_ligand {frag['name']} 无法解析"
        n_dummy = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 0)
        assert n_dummy == 1, f"e3_ligand {frag['name']} 应该有1个连接点，实际{n_dummy}"


def test_e3_ligand_molecular_weight_matches_known_crbn_ligand_scale():
    """
    来那度胺/泊马度胺是已上市药物，分子量是公开已知的(约259/273)——这里验证
    N-连接改造后的片段(减去被取代的那个H)分子量仍在同一量级，不是凭空编的结构。
    """
    libs = load_bifunctional_fragments(CONFIG_DIR)
    mw_by_name = {}
    for frag in libs["e3_ligands"]:
        mol = Chem.MolFromSmiles(frag["smiles"])
        mw_by_name[frag["name"]] = Descriptors.MolWt(mol)
    assert 250 < mw_by_name["lenalidomide_n_linked"] < 265
    assert 265 < mw_by_name["pomalidomide_n_linked"] < 280
