"""
三路线(route_a_shp2_sos1 / route_b_degrader / route_c_4th_gen_tki)共用同一套
core/ 计算引擎，区别只在于各自 routes/<route_id>/config/ 下的配置文件。

这个模块只做一件事：把 route_id 解析成配置目录路径。不做任何缓存/校验业务逻辑，
保持和其它 core/*.py 模块一样的"显式参数，不藏全局状态"风格——调用方自己决定
要用哪个route，不靠环境变量这类隐式开关。
"""
from __future__ import annotations

from pathlib import Path

ROUTES_DIR = Path(__file__).resolve().parent.parent / "routes"

DEFAULT_ROUTE = "route_c_4th_gen_tki"  # 历史默认值，保证所有现有调用方不传参数时行为不变


def route_config_dir(route_id: str = DEFAULT_ROUTE) -> Path:
    """返回 routes/<route_id>/config 目录，route_id不存在时明确报错而不是悄悄用错配置。"""
    d = ROUTES_DIR / route_id / "config"
    if not d.is_dir():
        available = sorted(p.name for p in ROUTES_DIR.iterdir() if p.is_dir()) if ROUTES_DIR.is_dir() else []
        raise FileNotFoundError(f"未知route或配置目录不存在: {d}（可用route: {available}）")
    return d


if __name__ == "__main__":
    print("=== 三条路线的配置目录 ===")
    for route_id in ("route_a_shp2_sos1", "route_b_degrader", "route_c_4th_gen_tki"):
        print(f"  {route_id} -> {route_config_dir(route_id)}")

    print("\n=== 未知route应该明确报错 ===")
    try:
        route_config_dir("route_z_does_not_exist")
    except FileNotFoundError as exc:
        print(f"  正确抛出: {exc}")
