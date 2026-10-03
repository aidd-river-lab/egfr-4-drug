#!/bin/bash
# 搭建一套能跑真实AutoDock Vina对接的独立环境(.venv310)，不碰原来的.venv(Python 3.9)。
#
# 为什么需要这个脚本，不能直接 pip install vina：
#   1. meeko/vina 需要 Python >= 3.10(主.venv是3.9.6装不了)
#   2. vina的setup.py硬编码只在 /usr/local/include 等几个固定路径找Boost，
#      Apple Silicon上Homebrew装在/opt/homebrew，不在那几个路径里——
#      用CONDA_PREFIX环境变量指向Homebrew的boost前缀骗过它的检测逻辑
#      (它的locate_boost()函数真的会去读这个环境变量，不是瞎编的trick)
#   3. vina 1.2.7的C++代码写死用 -std=c++11 编译，但较新的Boost(1.92+)的
#      <boost/math/tools/type_traits.hpp> 用到了C++14才有的std::make_signed_t
#      等类型别名，在严格C++11模式下编译会报错——下载vina源码后把这一个
#      编译选项从c++11改成c++17，问题就解决了(已经实测验证过)。
#   4. 还需要swig(生成Python绑定)，Homebrew没有默认装。
#
# 用法：bash scripts/setup_docking_env.sh
set -euo pipefail
cd "$(dirname "$0")/.."

echo "=== 1. Homebrew依赖 ==="
brew list python@3.10 >/dev/null 2>&1 || brew install python@3.10
brew list boost >/dev/null 2>&1 || brew install boost
brew list swig >/dev/null 2>&1 || brew install swig

echo "=== 2. 建 .venv310 ==="
/opt/homebrew/bin/python3.10 -m venv .venv310
.venv310/bin/pip install -q --upgrade pip
.venv310/bin/pip install -q rdkit pandas numpy pyyaml scipy gemmi meeko

echo "=== 3. 下载+打补丁+编译安装 vina(标准pip源没有macOS wheel，必须从源码编译) ==="
VINA_VERSION="1.2.7"
WORKDIR=$(mktemp -d)
curl -sL -o "$WORKDIR/vina.tar.gz" \
  "https://files.pythonhosted.org/packages/d2/2a/6746ef5e57b1c643e9fb24ad9e4fa520add7338736d50954e0fbc12ae52e/vina-${VINA_VERSION}.tar.gz"
tar xzf "$WORKDIR/vina.tar.gz" -C "$WORKDIR"
sed -i '' 's/"-std=c++11",/"-std=c++17",/' "$WORKDIR/vina-${VINA_VERSION}/setup.py"

export CONDA_DEFAULT_ENV=fake_env_for_boost_detection
export CONDA_PREFIX=$(brew --prefix boost)
.venv310/bin/pip install -q "$WORKDIR/vina-${VINA_VERSION}"
rm -rf "$WORKDIR"

echo "=== 4. 安装 PyRosetta(点突变建模用，core/structures.py::mutate_residue()) ==="
# pyrosetta-installer内部用subprocess调用裸`pip`(不是`sys.executable -m pip`)，
# 如果PATH上排在前面的pip指向别的Python版本(比如装了miniconda)，会把PyRosetta
# 装错环境、拿到不匹配的wheel报"not a supported wheel on this platform"——
# 显式把.venv310/bin塞到PATH最前面，强制它用这个venv自己的pip。
.venv310/bin/pip install -q pyrosetta-installer
PATH="$(pwd)/.venv310/bin:$PATH" .venv310/bin/python -c "
import pyrosetta_installer
pyrosetta_installer.install_pyrosetta()
"

echo "=== 5. pymysql(写DB用) ==="
.venv310/bin/pip install -q pymysql

echo "=== 6. 验证 ==="
.venv310/bin/python -c "
from vina import Vina
import meeko
import pyrosetta
pyrosetta.init('-mute all')
print('meeko + vina + pyrosetta 都能正常import')
"

echo ""
echo "完成。用法示例："
echo "  .venv310/bin/mk_prepare_receptor.py --read_pdb xxx.pdb -o receptor --allow_bad_res --default_altloc A -p"
echo "  .venv310/bin/mk_prepare_ligand.py -i xxx.sdf -o ligand.pdbqt"
echo "  .venv310/bin/python -c \"from core.docking import run_vina_docking; ...\""
echo "  .venv310/bin/python -c \"from core.structures import mutate_residue; ...\""
