"""core/route_config.py 的路径解析测试。"""
import pytest

from core.route_config import route_config_dir


@pytest.mark.parametrize("route_id", ["route_a_shp2_sos1", "route_b_degrader", "route_c_4th_gen_tki"])
def test_all_three_routes_resolve_to_existing_directory(route_id):
    d = route_config_dir(route_id)
    assert d.is_dir()
    assert d.name == "config"
    assert d.parent.name == route_id


def test_default_route_is_route_c():
    d = route_config_dir()
    assert d.parent.name == "route_c_4th_gen_tki"


def test_unknown_route_raises_with_available_list():
    with pytest.raises(FileNotFoundError) as exc_info:
        route_config_dir("route_z_does_not_exist")
    message = str(exc_info.value)
    assert "route_z_does_not_exist" in message
    assert "route_a_shp2_sos1" in message  # 错误信息要把可用route列出来，不能让人瞎猜
