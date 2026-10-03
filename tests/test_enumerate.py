"""
core/enumerate.py 的回归测试，重点锁死开发时真实踩过的那个molzip同位素不匹配bug：
裸 "*" 片段如果不强制改成骨架要求的同位素编号，molzip会静默留下断裂的多组分分子，
SanitizeMol不会报错——必须靠显式的连通性检查才能抓到。
"""
from rdkit import Chem

from core.enumerate import assemble, enumerate_from_scaffold


SCAFFOLD = "[1*]c1ccc(Nc2nccc(-c3ccc([2*])cc3[3*])n2)cc1"


def test_assemble_connects_all_fragments_into_single_component():
    smiles = assemble(
        SCAFFOLD,
        {1: "*S(C)(=O)=O", 2: "*C", 3: "*N1CCOCC1"},
    )
    assert smiles is not None
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None
    assert len(Chem.GetMolFrags(mol)) == 1


def test_assemble_leaves_no_dummy_atoms():
    smiles = assemble(
        SCAFFOLD,
        {1: "*S(C)(=O)=O", 2: "*C", 3: "*N1CCOCC1"},
    )
    mol = Chem.MolFromSmiles(smiles)
    assert not any(a.GetAtomicNum() == 0 for a in mol.GetAtoms())


def test_assemble_rejects_mismatched_isotope_as_broken_not_silent():
    """
    回归锁定：如果片段的哑原子同位素和骨架要求的不匹配(这里故意不传force_isotope，
    模拟旧bug)，assemble必须返回None，而不是悄悄吐出一个多组分SMILES。
    """
    core = Chem.MolFromSmiles(SCAFFOLD)
    frag = Chem.MolFromSmiles("*S(C)(=O)=O")  # isotope=0，骨架期望的是1/2/3
    combined = Chem.CombineMols(core, frag)
    params = Chem.MolzipParams()
    params.label = Chem.MolzipLabel.Isotope
    result = Chem.molzip(combined, params)
    Chem.SanitizeMol(result)  # 这一步不会报错，这正是当初这个bug容易被漏过的原因
    assert len(Chem.GetMolFrags(result)) > 1  # 证明确实断裂了，不是assemble()本身在骗人


def test_enumerate_from_scaffold_produces_fully_connected_unique_molecules():
    df = enumerate_from_scaffold("demo_aminopyrimidine_biphenyl")
    assert len(df) > 0
    assert df["smiles"].is_unique
    for smiles in df["smiles"]:
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None
        assert len(Chem.GetMolFrags(mol)) == 1
        assert not any(a.GetAtomicNum() == 0 for a in mol.GetAtoms())


def test_enumerate_count_matches_cartesian_product_of_region_sizes():
    df = enumerate_from_scaffold("demo_aminopyrimidine_biphenyl")
    # demo骨架: R1=entrance_797(6个), R2=hydrophobic_back_pocket(6个), R3=solvent_exposed(9个)
    assert len(df) == 6 * 6 * 9
