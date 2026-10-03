"""
config/structural_alerts.yaml 里每条规则的正/负对照验证。
这是 structural_alerts.yaml 顶部注释承诺的那份测试脚本——规则改动后必须跑一遍这里，
不能只靠眼看SMARTS语法"看起来对"就合入。
"""
from rdkit import Chem

from core.admet import run_structural_alerts

OSIMERTINIB = "COc1cc(N(C)CCN(C)C)c(NC(=O)C=C)cc1Nc1nccc(-c2cn(C)c3ccccc23)n1"
GEFITINIB = "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1"


def _hit_names(smiles, warhead_smarts=None):
    report = run_structural_alerts(smiles, warhead_smarts=warhead_smarts)
    return {h.rule_name for h in report.hits}


def test_osimertinib_no_hits_when_warhead_exempted():
    hits = _hit_names(OSIMERTINIB, warhead_smarts="[CX3]=[CX3][CX3]=[OX1]")
    assert hits == set()


def test_osimertinib_only_hits_michael_acceptor_when_not_exempted():
    hits = _hit_names(OSIMERTINIB)
    assert hits == {"michael_acceptor_offtarget"}


def test_gefitinib_aminoquinazoline_not_flagged_as_unmasked_aniline():
    mol = Chem.MolFromSmiles(GEFITINIB)
    assert mol is not None
    hits = _hit_names(GEFITINIB)
    assert "aniline_unmasked" not in hits


def test_plain_aniline_flagged():
    assert "aniline_unmasked" in _hit_names("Nc1ccccc1")


def test_acylated_aniline_not_flagged():
    assert "aniline_unmasked" not in _hit_names("CC(=O)Nc1ccccc1")


def test_para_halogenated_aniline_exempted():
    assert "aniline_unmasked" not in _hit_names("Nc1ccc(F)cc1")
    assert "aniline_unmasked" not in _hit_names("Nc1ccc(Cl)cc1")


def test_nitrobenzene_flagged_reject():
    hits = _hit_names("c1ccccc1[N+](=O)[O-]")
    assert "nitroaromatic" in hits


def test_aliphatic_nitro_not_flagged():
    assert "nitroaromatic" not in _hit_names("CC[N+](=O)[O-]")


def test_unsubstituted_furan_flagged():
    assert "furan_alpha_unsubstituted" in _hit_names("c1ccoc1")


def test_methylfuran_not_flagged():
    assert "furan_alpha_unsubstituted" not in _hit_names("Cc1ccco1")


def test_unsubstituted_thiophene_flagged():
    assert "thiophene_alpha_unsubstituted" in _hit_names("c1ccsc1")


def test_methylthiophene_not_flagged():
    assert "thiophene_alpha_unsubstituted" not in _hit_names("Cc1sccc1")


def test_free_pyrrole_flagged():
    assert "pyrrole_alpha_unsubstituted" in _hit_names("c1cc[nH]c1")


def test_n_methylpyrrole_not_flagged():
    assert "pyrrole_alpha_unsubstituted" not in _hit_names("Cn1cccc1")


def test_acrylamide_flagged_as_michael_acceptor():
    assert "michael_acceptor_offtarget" in _hit_names("C=CC(=O)N")


def test_para_phenylenediamine_flagged_reject():
    assert "diaminobenzene_1_4" in _hit_names("Nc1ccc(N)cc1")


def test_mono_acylated_para_phenylenediamine_not_flagged():
    assert "diaminobenzene_1_4" not in _hit_names("CC(=O)Nc1ccc(N)cc1")


def test_ortho_phenylenediamine_flagged_reject():
    assert "diaminobenzene_1_2" in _hit_names("Nc1ccccc1N")


def test_2_methylpyridine_flagged():
    assert "methylheterocycle_alpha_n" in _hit_names("Cc1ccccn1")


def test_3_methylpyridine_not_flagged():
    assert "methylheterocycle_alpha_n" not in _hit_names("Cc1cccnc1")


def test_unparseable_smiles_reports_parse_error_reject():
    report = run_structural_alerts("not_a_smiles(((")
    assert report.verdict == "REJECT"
    assert report.hits[0].rule_name == "parse_error"
