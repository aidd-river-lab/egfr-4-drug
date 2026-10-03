"""
core/enumerate.py 的 enumerate_bifunctional()/enumerate_bifunctional_library() 回归测试。
三组分(warhead-linker-E3)拓扑和固定骨架+R基团拓扑共用同一类molzip连通性风险，
这里用同样的验证强度覆盖：正常路径、断裂检测、端到端库枚举。
"""
from rdkit import Chem

from core.enumerate import enumerate_bifunctional, enumerate_bifunctional_library, load_bifunctional_fragments
from core.route_config import route_config_dir

WARHEAD = "*c1cc(C(F)(F)F)ccc1-c1ccnc(Nc2ccc(S(C)(=O)=O)cc2)n1"
LINKER = "*CCOCCOCC*"
E3_LIGAND = "*Nc1cccc2c1CN(C1CCC(=O)NC1=O)C2=O"


def test_enumerate_bifunctional_connects_all_three_fragments():
    smiles = enumerate_bifunctional(WARHEAD, LINKER, E3_LIGAND)
    assert smiles is not None
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None
    assert len(Chem.GetMolFrags(mol)) == 1


def test_enumerate_bifunctional_leaves_no_dummy_atoms():
    smiles = enumerate_bifunctional(WARHEAD, LINKER, E3_LIGAND)
    mol = Chem.MolFromSmiles(smiles)
    assert not any(a.GetAtomicNum() == 0 for a in mol.GetAtoms())


def test_enumerate_bifunctional_rejects_linker_with_wrong_dummy_count():
    """linker必须恰好2个连接点，给一个只有1个连接点的片段当linker应该诚实失败而不是产出断裂分子。"""
    single_dummy_fragment = "*CCCC"  # 只有一端有连接点
    result = enumerate_bifunctional(WARHEAD, single_dummy_fragment, E3_LIGAND)
    assert result is None


def test_load_bifunctional_fragments_from_route_b_config():
    libs = load_bifunctional_fragments(route_config_dir("route_b_degrader"))
    assert len(libs["warheads"]) >= 1
    assert len(libs["linkers"]) >= 1
    assert len(libs["e3_ligands"]) >= 1


def test_enumerate_bifunctional_library_produces_fully_connected_unique_molecules():
    df = enumerate_bifunctional_library(config_dir=route_config_dir("route_b_degrader"))
    assert len(df) > 0
    assert df["smiles"].is_unique
    for smiles in df["smiles"]:
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None
        assert len(Chem.GetMolFrags(mol)) == 1
        assert not any(a.GetAtomicNum() == 0 for a in mol.GetAtoms())


def test_enumerate_bifunctional_library_count_matches_cartesian_product():
    df = enumerate_bifunctional_library(config_dir=route_config_dir("route_b_degrader"))
    libs = load_bifunctional_fragments(route_config_dir("route_b_degrader"))
    expected = len(libs["warheads"]) * len(libs["linkers"]) * len(libs["e3_ligands"])
    assert len(df) == expected
