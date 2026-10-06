"""
core/structures.py 的测试。fetch_pdb/clean_chain是真实网络IO，不适合在每次pytest
里都跑(会打RCSB)，已经在本模块的__main__里验证过。这里只测mutate_residue()的
lazy import诚实失败路径——主.venv(Python 3.9)装不了PyRosetta，调用时应该诚实
返回ok=False，不是抛异常崩溃，也不是假装成功。
"""
from core.structures import mutate_residue, run_protein_equilibration_md


def test_mutate_residue_honestly_fails_without_pyrosetta():
    """
    本仓库的主.venv没有装PyRosetta(真实调用需要.venv310，见
    scripts/setup_docking_env.sh)，所以这个测试在CI/主venv下预期ok=False。
    如果哪天主venv也装了PyRosetta，这个测试会跑到真实分支——那也是好事，
    不需要改这个测试本身(mutate_residue用的是文件/链/残基号这种通用参数，
    不是"假设没装"这种脆弱的mock)。
    """
    result = mutate_residue(
        structure_path="/nonexistent/path.pdb", chain="A", resnum=790,
        new_aa_one_letter="T", out_path="/tmp/should_not_be_created.pdb",
    )
    if not result["ok"]:
        assert "PyRosetta" in result["reason"]
    else:
        # 如果运行环境确实装了PyRosetta，说明上面的文件路径找不到应该在pyrosetta内部报错，
        # 不应该意外成功
        raise AssertionError("不存在的文件路径不应该成功返回ok=True")


def test_run_protein_equilibration_md_honestly_fails_without_openmm():
    """
    主.venv(Python 3.9)没有装openmm/pdbfixer(真实调用需要.venv310，2026-10起
    已经真实装好并验证跑通，见doc/02-structure-ensemble.md)。同样的"不是假设
    没装"原则：这个测试检查的是lazy import失败时的诚实返回，不是mock掉openmm。
    """
    result = run_protein_equilibration_md(
        structure_path="/nonexistent/path.pdb",
        out_trajectory_path="/tmp/should_not_be_created.dcd",
    )
    if not result["ok"]:
        assert "openmm" in result["reason"]
    else:
        raise AssertionError("不存在的文件路径不应该成功返回ok=True")
