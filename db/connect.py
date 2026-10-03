"""
MySQL连接helper，从项目根目录的 .env 读取凭证(不读环境变量以外的任何硬编码值，
真实host/用户名/密码绝不出现在代码里——这个文件本身应该被git跟踪，.env不应该)。

用法：
    from db.connect import get_connection
    conn = get_connection()               # 连到 .env 里 DB_NAME 指定的库
    conn = get_connection(database=None)  # 不指定库，用于建库这类操作
"""
from __future__ import annotations

from pathlib import Path

import pymysql

ENV_PATH = Path(__file__).resolve().parent.parent / ".env"


def _load_env() -> dict:
    if not ENV_PATH.is_file():
        raise FileNotFoundError(
            f"{ENV_PATH} 不存在。复制 .env.example 为 .env 并填入真实的数据库凭证(不要提交到git)。"
        )
    env = {}
    with open(ENV_PATH, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#") and "=" in line:
                key, value = line.split("=", 1)
                env[key.strip()] = value.strip()
    return env


def get_connection(database: str | None = "__default__", autocommit: bool = True):
    """
    database="__default__"(默认) 用 .env 里的 DB_NAME；传 None 则不指定库(建库场景)；
    也可以显式传别的库名。
    """
    env = _load_env()
    db_name = env.get("DB_NAME") if database == "__default__" else database
    return pymysql.connect(
        host=env["DB_HOST"],
        port=int(env.get("DB_PORT", 3306)),
        user=env["DB_USER"],
        password=env["DB_PASSWORD"],
        database=db_name,
        autocommit=autocommit,
        connect_timeout=10,
        charset="utf8mb4",
    )


if __name__ == "__main__":
    conn = get_connection()
    with conn.cursor() as cur:
        cur.execute("SELECT VERSION()")
        print("MySQL version:", cur.fetchone()[0])
        cur.execute("SHOW TABLES")
        tables = [r[0] for r in cur.fetchall()]
        print(f"{len(tables)} 张表:", sorted(tables))
    conn.close()
