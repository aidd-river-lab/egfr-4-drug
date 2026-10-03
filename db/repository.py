"""
把 workflows/*.py 算出来的DataFrame/dict结果写进MySQL——这一层只做数据搬运，
不做任何业务计算(业务逻辑都在core/里)，也不决定"该不该算"，只决定"算完的东西
怎么落库"。

幂等策略(重跑同一轮workflow会发生什么)：
  - compounds/designed_molecules/admet_descriptors/synthesis_records：这些是
    "当前状态"类的表(每个compound_id一行)，用 INSERT...ON DUPLICATE KEY UPDATE，
    重跑会用最新计算结果覆盖旧值，不会堆积重复行。
  - structural_alert_hits：先删除涉及到的compound_id的旧记录，再插入新的——
    同样是"当前状态"语义(重跑后警示状态应该反映最新规则/最新结构)，不是日志表。
  - standardize_quarantine：纯追加，没有唯一键——隔离记录本质是审计日志，
    每次失败都值得留痕，不需要去重。
  - decision_rounds/decision_batch_members：每次调用生成一个新的round_id
    (路线名+时间戳)，代表"这是哪一轮DMTA决策"，天然不会和之前的round冲突。

每个 upsert_* / insert_* 函数内部调用一个同名前缀 _*_rows() 的纯函数做
"DataFrame -> SQL参数元组列表"的转换，这部分不需要真实数据库连接就能单测
(见 tests/test_repository.py)，和项目里其它"纯逻辑先独立验证"的模式一致。
"""
from __future__ import annotations

import json
from datetime import date, datetime

import pandas as pd


# ------------------------------------------------------------
# 环节1：compounds / standardize_quarantine
# ------------------------------------------------------------
def _compound_rows(df: pd.DataFrame, route_id: str) -> list[tuple]:
    rows = []
    for _, r in df.iterrows():
        rows.append(
            (
                r["compound_id"],
                route_id,
                r["standardized_smiles"],
                r["inchikey"],
                r["inchikey14"],
                r.get("source", "enumeration"),
                r.get("scaffold_id"),
            )
        )
    return rows


class CrossRouteCollisionError(Exception):
    """compound_id在compounds表里已经存在，但属于另一条路线——绝不能静默覆盖。"""


def _check_no_cross_route_collision(conn, compound_ids: list[str], route_id: str) -> None:
    """
    真实事故复盘：曾经因为compound_id生成逻辑按scaffold_id截断(而不是按route_id)，
    导致路线A和路线C生成了同一批compound_id，upsert时(PRIMARY KEY只有compound_id，
    没有route_id)互相覆盖了对方的smiles_canonical——route_id字段显示的是先插入那条
    路线的，但分子数据却是后插入那条路线的，是一次真实的跨路线数据污染。
    根因已经在 workflows/run_discovery_round.py 里修了(compound_id前缀改成按route_id生成，
    不同路线永远不会撞ID)，但这里额外加一道防线：upsert前检查，如果某个compound_id
    已经存在但属于别的route_id，直接报错而不是让ON DUPLICATE KEY UPDATE静默覆盖。
    """
    if not compound_ids:
        return
    placeholders = ",".join(["%s"] * len(compound_ids))
    with conn.cursor() as cur:
        cur.execute(
            f"SELECT compound_id, route_id FROM compounds WHERE compound_id IN ({placeholders}) AND route_id != %s",
            (*compound_ids, route_id),
        )
        conflicts = cur.fetchall()
    if conflicts:
        detail = ", ".join(f"{cid}(属于{existing_route})" for cid, existing_route in conflicts[:5])
        more = f" 等共{len(conflicts)}个" if len(conflicts) > 5 else ""
        raise CrossRouteCollisionError(
            f"拒绝写入：以下compound_id已经存在但属于别的路线，不能被route_id={route_id}的数据覆盖——"
            f"{detail}{more}。这说明compound_id生成逻辑在不同路线之间产生了冲突，去检查ID生成规则，"
            f"不要在这里强行覆盖。"
        )


def upsert_compounds(conn, df: pd.DataFrame, route_id: str) -> int:
    rows = _compound_rows(df, route_id)
    if not rows:
        return 0
    _check_no_cross_route_collision(conn, [r[0] for r in rows], route_id)
    sql = """
        INSERT INTO compounds (compound_id, route_id, smiles_canonical, inchikey, inchikey_skeleton, source, scaffold_id)
        VALUES (%s, %s, %s, %s, %s, %s, %s)
        ON DUPLICATE KEY UPDATE
            smiles_canonical = VALUES(smiles_canonical),
            inchikey_skeleton = VALUES(inchikey_skeleton),
            source = VALUES(source),
            scaffold_id = VALUES(scaffold_id)
    """
    with conn.cursor() as cur:
        cur.executemany(sql, rows)
    return len(rows)


def _quarantine_rows(df: pd.DataFrame) -> list[tuple]:
    rows = []
    for _, r in df.iterrows():
        reason = f"{r.get('fail_stage', '?')}: {r.get('fail_reason', '?')}"
        rows.append((r["raw_smiles"], reason))
    return rows


def insert_quarantine(conn, df: pd.DataFrame) -> int:
    if df.empty:
        return 0
    rows = _quarantine_rows(df)
    with conn.cursor() as cur:
        cur.executemany("INSERT INTO standardize_quarantine (raw_smiles, reason) VALUES (%s, %s)", rows)
    return len(rows)


# ------------------------------------------------------------
# 环节3：designed_molecules
# ------------------------------------------------------------
def _designed_molecule_rows(
    df: pd.DataFrame, scaffold_id: str, enumeration_batch: str, r_group_cols: list[str]
) -> list[tuple]:
    rows = []
    for _, r in df.iterrows():
        assignment = {
            col[:-len("_name")]: r[col] for col in r_group_cols if col in r.index and pd.notna(r[col])
        }
        rows.append((r["compound_id"], scaffold_id, json.dumps(assignment, ensure_ascii=False), enumeration_batch))
    return rows


def upsert_designed_molecules(
    conn, df: pd.DataFrame, scaffold_id: str, enumeration_batch: str, r_group_cols: list[str]
) -> int:
    rows = _designed_molecule_rows(df, scaffold_id, enumeration_batch, r_group_cols)
    if not rows:
        return 0
    sql = """
        INSERT INTO designed_molecules (compound_id, scaffold_id, r_group_assignment, enumeration_batch)
        VALUES (%s, %s, %s, %s)
        ON DUPLICATE KEY UPDATE
            r_group_assignment = VALUES(r_group_assignment),
            enumeration_batch = VALUES(enumeration_batch)
    """
    with conn.cursor() as cur:
        cur.executemany(sql, rows)
    return len(rows)


# ------------------------------------------------------------
# 环节6：admet_descriptors / structural_alert_hits
# ------------------------------------------------------------
def _bool_or_none(value):
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    return int(bool(value))


def _admet_rows(df: pd.DataFrame) -> list[tuple]:
    rows = []
    for _, r in df.iterrows():
        rows.append(
            (
                r["compound_id"],
                r.get("mw"),
                r.get("clogp"),
                r.get("clogd_approx"),
                r.get("tpsa"),
                r.get("hbd"),
                r.get("hba"),
                r.get("fsp3"),
                _bool_or_none(r.get("has_basic_aliphatic_amine")),
                r.get("pka_proxy"),
                r.get("cns_mpo"),
                _bool_or_none(r.get("cns_mpo_pass_4_0")),
                r.get("herg_risk_proxy"),
            )
        )
    return rows


def upsert_admet_descriptors(conn, df: pd.DataFrame) -> int:
    rows = _admet_rows(df)
    if not rows:
        return 0
    sql = """
        INSERT INTO admet_descriptors
            (compound_id, mw, clogp, clogd_approx, tpsa, hbd, hba, fsp3,
             has_basic_aliphatic_amine, pka_proxy, cns_mpo, cns_mpo_pass_4_0, herg_risk_proxy)
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
        ON DUPLICATE KEY UPDATE
            mw = VALUES(mw), clogp = VALUES(clogp), clogd_approx = VALUES(clogd_approx), tpsa = VALUES(tpsa),
            hbd = VALUES(hbd), hba = VALUES(hba), fsp3 = VALUES(fsp3),
            has_basic_aliphatic_amine = VALUES(has_basic_aliphatic_amine), pka_proxy = VALUES(pka_proxy),
            cns_mpo = VALUES(cns_mpo), cns_mpo_pass_4_0 = VALUES(cns_mpo_pass_4_0),
            herg_risk_proxy = VALUES(herg_risk_proxy)
    """
    with conn.cursor() as cur:
        cur.executemany(sql, rows)
    return len(rows)


def replace_alert_hits(conn, compound_ids: list[str], hit_rows: list[tuple]) -> int:
    """
    hit_rows: (compound_id, rule_name, severity, action, advice) 元组列表。
    先删除这些compound_id已有的命中记录，再插入新的，保持"当前警示状态"语义
    (见模块docstring)，不是在tests/test_repository.py里的_*_rows()那类纯转换函数——
    这个函数天生需要先DELETE再INSERT两步，拆成纯函数意义不大，直接在这里做。
    """
    with conn.cursor() as cur:
        if compound_ids:
            placeholders = ",".join(["%s"] * len(compound_ids))
            cur.execute(f"DELETE FROM structural_alert_hits WHERE compound_id IN ({placeholders})", compound_ids)
        if hit_rows:
            cur.executemany(
                "INSERT INTO structural_alert_hits (compound_id, rule_name, severity, action, advice) "
                "VALUES (%s, %s, %s, %s, %s)",
                hit_rows,
            )
    return len(hit_rows)


# ------------------------------------------------------------
# 环节7：synthesis_records(目前只有sa_score，其它字段诚实留空，见core/synthesis.py)
# ------------------------------------------------------------
def _synthesis_rows(df: pd.DataFrame) -> list[tuple]:
    return [(r["compound_id"], r.get("sa_score")) for _, r in df.iterrows()]


def upsert_synthesis_records(conn, df: pd.DataFrame) -> int:
    rows = _synthesis_rows(df)
    if not rows:
        return 0
    sql = """
        INSERT INTO synthesis_records (compound_id, sa_score) VALUES (%s, %s)
        ON DUPLICATE KEY UPDATE sa_score = VALUES(sa_score)
    """
    with conn.cursor() as cur:
        cur.executemany(sql, rows)
    return len(rows)


# ------------------------------------------------------------
# 环节8：decision_rounds / decision_batch_members(路线B目前没有决策阶段，不调用这两个)
# ------------------------------------------------------------
def _decision_batch_member_rows(round_id: str, df: pd.DataFrame) -> list[tuple]:
    rows = []
    for _, r in df.iterrows():
        rows.append(
            (
                round_id,
                r["compound_id"],
                r["batch_role"],
                _bool_or_none(r.get("is_pareto_optimal")),
                r.get("desirability"),
                None,  # hypothesis_label：本仓库的批次选择没有药化指定的假设文本，诚实留空
            )
        )
    return rows


def insert_decision_round(conn, route_id: str, batch_size: int, quota: dict, notes: str = "") -> str:
    round_id = f"{route_id[:20]}-{datetime.now():%Y%m%d%H%M%S}"
    with conn.cursor() as cur:
        cur.execute(
            "INSERT INTO decision_rounds (round_id, round_date, batch_size, quota_json, notes) VALUES (%s,%s,%s,%s,%s)",
            (round_id, date.today().isoformat(), batch_size, json.dumps(quota, ensure_ascii=False), notes),
        )
    return round_id


def insert_decision_batch_members(conn, round_id: str, df: pd.DataFrame) -> int:
    rows = _decision_batch_member_rows(round_id, df)
    if not rows:
        return 0
    sql = """
        INSERT INTO decision_batch_members
            (round_id, compound_id, batch_role, is_pareto_optimal, desirability_score, hypothesis_label)
        VALUES (%s, %s, %s, %s, %s, %s)
    """
    with conn.cursor() as cur:
        cur.executemany(sql, rows)
    return len(rows)


# ------------------------------------------------------------
# 顶层编排：一个workflow跑完之后，一次调用把这一轮能落库的东西全部落库
# ------------------------------------------------------------
def persist_discovery_round(
    conn,
    route_id: str,
    scaffold_id: str,
    enumeration_batch: str,
    r_group_cols: list[str],
    curated_df: pd.DataFrame,
    quarantine_df: pd.DataFrame,
    all_candidates_df: pd.DataFrame,
    alert_hit_rows: list[tuple],
    batch_combined_df: pd.DataFrame | None = None,
    batch_size: int | None = None,
    quota: dict | None = None,
) -> dict:
    """
    curated_df/all_candidates_df 是标准化之后、REJECT硬门过滤**之前**的全量候选——
    被REJECT掉的分子也有真实计算出来的ADMET/警示/SA score，这些数据本身是有价值的
    (比如以后做"规则有没有误杀"的复盘)，不应该因为没进最终批次就不落库。

    batch_combined_df 不传时(目前路线B就是这种情况，决策阶段还没对接真实Kd数据)，
    不写 decision_rounds/decision_batch_members，不是遗漏，是没有东西可写。
    """
    summary = {
        "compounds": upsert_compounds(conn, curated_df, route_id),
        "quarantine": insert_quarantine(conn, quarantine_df),
        "designed_molecules": upsert_designed_molecules(conn, curated_df, scaffold_id, enumeration_batch, r_group_cols),
        "admet_descriptors": upsert_admet_descriptors(conn, all_candidates_df),
        "structural_alert_hits": replace_alert_hits(conn, all_candidates_df["compound_id"].tolist(), alert_hit_rows),
        "synthesis_records": upsert_synthesis_records(conn, all_candidates_df),
    }
    if batch_combined_df is not None:
        round_id = insert_decision_round(conn, route_id, batch_size, quota or {})
        summary["decision_round_id"] = round_id
        summary["decision_batch_members"] = insert_decision_batch_members(conn, round_id, batch_combined_df)
    return summary
