"""
core/synthesis.py 的回归测试——之前这个模块完全没有pytest覆盖(只在__main__里
演示过)，这里补齐。compute_sa_score是纯RDKit，check_fto是纯文件读取+规则判断，
都不需要外部工具，可以在主.venv下完整验证。run_retrosynthesis()需要真实
AiZynthFinder+约1.2GB的预训练模型(见该函数docstring)，这里只测主.venv下的
诚实失败路径；真实执行结果(奥希替尼4步合成，4个precursor全部在ZINC库存里)
记录见doc/07-synthesis.md。
"""
from core.synthesis import check_fto, compute_sa_score, run_retrosynthesis, run_retrosynthesis_stub

OSIMERTINIB = "COc1cc(N(C)CCN(C)C)c(NC(=O)C=C)cc1Nc1nccc(-c2cn(C)c3ccccc23)n1"


def test_compute_sa_score_on_real_osimertinib():
    result = compute_sa_score(OSIMERTINIB)
    assert result["ok"] is True
    assert 1 <= result["sa_score"] <= 10


def test_compute_sa_score_honestly_fails_on_invalid_smiles():
    result = compute_sa_score("not a valid smiles!!!")
    assert result["ok"] is False


def test_compute_sa_score_accepts_rdkit_mol_directly():
    from rdkit import Chem

    mol = Chem.MolFromSmiles("c1ccccc1")
    result = compute_sa_score(mol)
    assert result["ok"] is True


def test_run_retrosynthesis_stub_is_honest_about_not_running():
    result = run_retrosynthesis_stub(OSIMERTINIB)
    assert result["ok"] is False
    assert result["route_found"] is None


def test_run_retrosynthesis_honestly_fails_without_aizynthfinder_or_config():
    """
    主.venv(Python 3.9)没有装aizynthfinder(真实调用需要.venv310 + 单独下载的
    ~1.2GB模型数据，2026-10起已经真实装好并验证跑通，见doc/07-synthesis.md)。
    这里传一个必然不存在的config路径，预期诚实返回ok=False——不管是因为
    aizynthfinder本身没装，还是因为配置文件找不到，都不应该抛异常或假装成功。
    """
    result = run_retrosynthesis(OSIMERTINIB, aizynth_config_path="/nonexistent/config.yml")
    assert result["ok"] is False


def test_check_fto_returns_unknown_with_placeholder_patent_data():
    """三条路线的patent_landscape.yaml目前都是status=placeholder_not_real，
    check_fto必须诚实返回unknown，不能伪装成真实的FTO结论。"""
    result = check_fto(OSIMERTINIB)
    assert result["ok"] is True
    assert result["fto_status"] == "unknown_no_real_patent_data"
