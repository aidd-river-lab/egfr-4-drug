"""
db/repository.py 里"DataFrame -> SQL参数元组"这部分纯转换逻辑的单测，不需要真实
数据库连接(upsert_*/insert_*那层薄薄的cursor.executemany()包装不在这里测，
需要真实MySQL连接的验证见本次人工跑通run_discovery_round(persist=True)的记录，
不适合写进自动化测试——CI环境不会有这个仓库.env里的真实远程凭证)。
"""
import json
from unittest.mock import MagicMock

import pandas as pd
import pytest

from db.repository import (
    CrossRouteCollisionError,
    _admet_rows,
    _bool_or_none,
    _check_no_cross_route_collision,
    _compound_rows,
    _decision_batch_member_rows,
    _designed_molecule_rows,
    _quarantine_rows,
    _synthesis_rows,
)


def test_compound_rows_shape_and_values():
    df = pd.DataFrame(
        [{"compound_id": "C1", "standardized_smiles": "CCO", "inchikey": "KEY1", "inchikey14": "KEY1SHORT"}]
    )
    rows = _compound_rows(df, route_id="route_c_4th_gen_tki")
    assert rows == [("C1", "route_c_4th_gen_tki", "CCO", "KEY1", "KEY1SHORT", "enumeration", None)]


def test_compound_rows_respects_explicit_source_and_scaffold():
    df = pd.DataFrame(
        [
            {
                "compound_id": "C1",
                "standardized_smiles": "CCO",
                "inchikey": "KEY1",
                "inchikey14": "KEY1SHORT",
                "source": "wet_lab",
                "scaffold_id": "demo_scaffold",
            }
        ]
    )
    rows = _compound_rows(df, route_id="route_a_shp2_sos1")
    assert rows[0] == ("C1", "route_a_shp2_sos1", "CCO", "KEY1", "KEY1SHORT", "wet_lab", "demo_scaffold")


def test_quarantine_rows_combines_stage_and_reason():
    df = pd.DataFrame([{"raw_smiles": "bad(((", "fail_stage": "parse", "fail_reason": "invalid syntax"}])
    rows = _quarantine_rows(df)
    assert rows == [("bad(((", "parse: invalid syntax")]


def test_designed_molecule_rows_packs_r_groups_into_json():
    df = pd.DataFrame([{"compound_id": "C1", "R1_name": "methylsulfone", "R2_name": "methyl", "R3_name": None}])
    rows = _designed_molecule_rows(df, scaffold_id="demo_scaffold", enumeration_batch="batch1", r_group_cols=["R1_name", "R2_name", "R3_name"])
    assert len(rows) == 1
    compound_id, scaffold_id, assignment_json, batch = rows[0]
    assignment = json.loads(assignment_json)
    assert assignment == {"R1": "methylsulfone", "R2": "methyl"}  # R3是None，应该被跳过，不是写成null


def test_bool_or_none_handles_nan_true_false():
    assert _bool_or_none(None) is None
    assert _bool_or_none(float("nan")) is None
    assert _bool_or_none(True) == 1
    assert _bool_or_none(False) == 0
    assert _bool_or_none(1) == 1


def test_admet_rows_shape():
    df = pd.DataFrame(
        [
            {
                "compound_id": "C1",
                "mw": 300.5,
                "clogp": 2.1,
                "has_basic_aliphatic_amine": True,
                "cns_mpo_pass_4_0": False,
            }
        ]
    )
    rows = _admet_rows(df)
    assert len(rows) == 1
    assert rows[0][0] == "C1"
    assert rows[0][1] == 300.5
    assert rows[0][8] == 1  # has_basic_aliphatic_amine -> 1
    assert rows[0][11] == 0  # cns_mpo_pass_4_0 -> 0


def test_synthesis_rows_shape():
    df = pd.DataFrame([{"compound_id": "C1", "sa_score": 2.5}])
    assert _synthesis_rows(df) == [("C1", 2.5)]


def test_decision_batch_member_rows_shape():
    df = pd.DataFrame(
        [{"compound_id": "C1", "batch_role": "exploit", "is_pareto_optimal": True, "desirability": 0.9}]
    )
    rows = _decision_batch_member_rows("round-1", df)
    assert rows == [("round-1", "C1", "exploit", 1, 0.9, None)]


def _mock_conn(existing_rows):
    """existing_rows: 模拟SELECT查到的(compound_id, route_id)冲突记录。"""
    conn = MagicMock()
    cursor = MagicMock()
    cursor.fetchall.return_value = existing_rows
    conn.cursor.return_value.__enter__.return_value = cursor
    return conn


def test_cross_route_collision_detected():
    """
    回归锁定一次真实的数据污染事故：compound_id生成逻辑曾经按scaffold_id截断，
    导致路线A和路线C生成了同一批compound_id，upsert时互相覆盖了对方的分子数据。
    这个检查函数必须在发现compound_id已属于别的route_id时报错，而不是放行。
    """
    conn = _mock_conn([("RTC-0000", "route_c_4th_gen_tki")])
    with pytest.raises(CrossRouteCollisionError, match="RTC-0000"):
        _check_no_cross_route_collision(conn, ["RTC-0000"], route_id="route_a_shp2_sos1")


def test_no_collision_when_same_route_reinserts():
    """同一条路线重复插入自己的compound_id(重跑workflow)不应该被当成冲突——
    _check_no_cross_route_collision 的SQL本身用 route_id != %s 过滤掉了这种情况，
    这里用一个空的existing_rows模拟查询结果(因为SQL层已经排除了同route的记录)。
    """
    conn = _mock_conn([])
    _check_no_cross_route_collision(conn, ["RTC-0000"], route_id="route_c_4th_gen_tki")  # 不应该抛异常


def test_no_collision_check_when_no_compound_ids():
    conn = _mock_conn([])
    _check_no_cross_route_collision(conn, [], route_id="route_c_4th_gen_tki")
    conn.cursor.assert_not_called()  # 空列表时不应该多发一次查询
