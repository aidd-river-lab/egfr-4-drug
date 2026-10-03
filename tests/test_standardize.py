"""
core/standardize.py 回归测试，重点锁死一个端到端集成才暴露出来的bug：
全部成功(或全部失败)时，DataFrame里压根不存在fail_stage/fail_reason列，
按列名取子集会抛KeyError——单独跑每个函数的demo不会触发，必须混合批量数据才会暴露。
"""
from core.standardize import standardize_batch

OSIMERTINIB = "COc1cc(N(C)CCN(C)C)c(NC(=O)C=C)cc1Nc1nccc(-c2cn(C)c3ccccc23)n1"


def test_all_success_batch_does_not_crash_on_empty_quarantine():
    result = standardize_batch([OSIMERTINIB, "CCO"])
    assert len(result.curated) == 2
    assert len(result.quarantine) == 0


def test_all_failure_batch_does_not_crash_on_empty_curated():
    result = standardize_batch(["not-a-smiles(((", "also-not-one((("])
    assert len(result.curated) == 0
    assert len(result.quarantine) == 2


def test_mixed_batch_splits_correctly():
    result = standardize_batch([OSIMERTINIB, "not-a-smiles((("])
    assert len(result.curated) == 1
    assert len(result.quarantine) == 1


def test_salt_stripped_and_dedup_key_present():
    result = standardize_batch(["c1ccccc1.[Na+].[Cl-]"])
    assert len(result.curated) == 1
    assert "inchikey" in result.curated.columns
    assert result.curated.loc[0, "inchikey"] is not None


def test_osimertinib_inchikey_matches_known_value():
    result = standardize_batch([OSIMERTINIB])
    assert result.curated.loc[0, "inchikey"] == "DUYJMQONPNNFPI-UHFFFAOYSA-N"
