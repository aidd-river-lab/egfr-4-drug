"""
core/covalent.py 的回归测试。evaluate_attack_geometry()已有的测试(合成
CovalentDockingResult数据)不重复；这里专门测新增的
extract_michael_acceptor_geometry()——从真实PDBQT格式的3D坐标里提取弹头
几何，纯几何计算，不需要PyRosetta/Vina，用手写的、坐标已知的合成PDBQT
文件验证距离/角度算对了没有。真实奥希替尼+PDB 4ZAU跑出来的结果见
doc/15-route-specific-gaps-and-next-steps.md。
"""
from core.covalent import extract_michael_acceptor_geometry


def _pdbqt_line(serial, atom_name, resname, resnum, x, y, z, atype):
    """按core/covalent.py::_parse_pdbqt_atoms()解析用的固定列位置，构造一行合法PDBQT文本。"""
    line = list(" " * 80)
    line[0:4] = list("ATOM")
    line[6:11] = list(f"{serial:>5}")
    line[12:16] = list(f"{atom_name:<4}")
    line[17:20] = list(f"{resname:<3}")
    line[22:26] = list(f"{resnum:>4}")
    line[30:38] = list(f"{x:8.3f}")
    line[38:46] = list(f"{y:8.3f}")
    line[46:54] = list(f"{z:8.3f}")
    line[77:79] = list(f"{atype:<2}")
    return "".join(line) + "\n"


def _write_ligand_with_acrylamide(path, extra_degree1_pairs=()):
    """构造一个含唯一丙烯酰胺端位烯烃碳(CB-CA键长1.321Å，C=C量级)的合成配体。
    extra_degree1_pairs: 额外的(bond_length, element)对，用来测试"长单键不该被误判成双键"
    这类边界情况，不破坏CB的唯一性。"""
    lines = [
        _pdbqt_line(1, "CB", "LIG", 1, 0.0, 0.0, 0.0, "C"),  # 端位烯烃碳，degree=1，和CA键长1.321
        _pdbqt_line(2, "CA", "LIG", 1, 1.321, 0.0, 0.0, "C"),  # CA还连着CO，degree=2，不会被误判
        _pdbqt_line(3, "CO", "LIG", 1, 2.600, 0.0, 0.0, "C"),  # 羰基碳
        # 羰基氧，让CO的degree变成2(真实丙烯酰胺里CO还连着=O和N，这里只加O就够让它
        # 不再是degree-1——漏了这个原子是本测试最初版本的bug，CO会被误判成第二个候选)
        _pdbqt_line(4, "OX", "LIG", 1, 2.600, 1.230, 0.0, "O"),
    ]
    serial = 5
    x_offset = 10.0
    for bond_length, element in extra_degree1_pairs:
        lines.append(_pdbqt_line(serial, f"X{serial}", "LIG", 1, x_offset, 0.0, 0.0, "C"))
        lines.append(
            _pdbqt_line(serial + 1, f"Y{serial}", "LIG", 1, x_offset + bond_length, 0.0, 0.0, element)
        )
        serial += 2
        x_offset += 10.0
    with open(path, "w") as f:
        f.writelines(lines)


def _write_receptor_with_nucleophile(path, x, y, z, resname="CYS", resnum=797, atom_name="SG"):
    with open(path, "w") as f:
        f.write(_pdbqt_line(1, atom_name, resname, resnum, x, y, z, "SA"))


def test_extract_geometry_known_distance_and_angle(tmp_path):
    """手算验证：Cβ在原点，Cα在(1.321,0,0)，亲核原子在(0,3,4)——
    距离应该是5.0(3-4-5直角三角形)，v1=(0,3,4)和v2=(1.321,0,0)正交，角度应该是90°。"""
    ligand = tmp_path / "ligand.pdbqt"
    receptor = tmp_path / "receptor.pdbqt"
    _write_ligand_with_acrylamide(ligand)
    _write_receptor_with_nucleophile(receptor, 0.0, 3.0, 4.0)

    result = extract_michael_acceptor_geometry(ligand, receptor, "CYS", 797, "SG")
    assert result["ok"] is True
    assert result["beta_carbon_to_nucleophile_distance_angstrom"] == 5.0
    assert result["attack_angle_degree"] == 90.0
    assert result["beta_alpha_bond_length_angstrom"] == 1.321


def test_extract_geometry_ignores_long_single_bonds(tmp_path):
    """额外加一对键长1.54Å(典型C-C单键)的碳——不应该被误判成端位烯烃碳，
    CB依然应该是唯一被识别出来的那个。"""
    ligand = tmp_path / "ligand.pdbqt"
    receptor = tmp_path / "receptor.pdbqt"
    _write_ligand_with_acrylamide(ligand, extra_degree1_pairs=[(1.54, "C")])
    _write_receptor_with_nucleophile(receptor, 0.0, 3.0, 4.0)

    result = extract_michael_acceptor_geometry(ligand, receptor, "CYS", 797, "SG")
    assert result["ok"] is True
    assert result["beta_carbon_to_nucleophile_distance_angstrom"] == 5.0


def test_extract_geometry_honestly_fails_when_ambiguous(tmp_path):
    """两组都符合"degree-1碳+碳邻居+键长<1.40"的候选——无法唯一确定，诚实返回ok=False。"""
    ligand = tmp_path / "ligand.pdbqt"
    receptor = tmp_path / "receptor.pdbqt"
    _write_ligand_with_acrylamide(ligand, extra_degree1_pairs=[(1.33, "C")])
    _write_receptor_with_nucleophile(receptor, 0.0, 3.0, 4.0)

    result = extract_michael_acceptor_geometry(ligand, receptor, "CYS", 797, "SG")
    assert result["ok"] is False


def test_extract_geometry_honestly_fails_without_nucleophile(tmp_path):
    """受体里根本没有指定的亲核原子(比如C797S突变体系，SG不存在)——诚实失败。"""
    ligand = tmp_path / "ligand.pdbqt"
    receptor = tmp_path / "receptor.pdbqt"
    _write_ligand_with_acrylamide(ligand)
    _write_receptor_with_nucleophile(receptor, 0.0, 3.0, 4.0, resname="SER", atom_name="OG")  # C797S，没有SG

    result = extract_michael_acceptor_geometry(ligand, receptor, "CYS", 797, "SG")
    assert result["ok"] is False
