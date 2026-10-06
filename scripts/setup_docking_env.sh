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

echo "=== 7. OpenMM + PDBFixer(纯pip装，蛋白短程MD平衡用，core/structures.py::run_protein_equilibration_md()) ==="
.venv310/bin/pip install -q openmm pdbfixer

echo "=== 7b. AiZynthFinder(纯pip装，真实逆合成路线搜索，core/synthesis.py::run_retrosynthesis()) ==="
.venv310/bin/pip install -q aizynthfinder
echo "    还需要下载真实预训练模型+ZINC建块库(约1.2GB，不随代码库提交)："
echo "    .venv310/bin/download_public_data validation/aizynth_data"

echo ""
echo "完成。用法示例："
echo "  .venv310/bin/mk_prepare_receptor.py --read_pdb xxx.pdb -o receptor --allow_bad_res --default_altloc A -p"
echo "  .venv310/bin/mk_prepare_ligand.py -i xxx.sdf -o ligand.pdbqt"
echo "  .venv310/bin/python -c \"from core.docking import run_vina_docking; ...\""
echo "  .venv310/bin/python -c \"from core.structures import mutate_residue; ...\""
echo "  .venv310/bin/python -c \"from core.structures import run_protein_equilibration_md; ...\""
echo ""
echo "注意：配体结合复合物的MD(core/md_stability.py::run_protein_ligand_complex_md())"
echo "需要openff-toolkit，它需要Python>=3.11或者conda解决一套复杂依赖，pip在.venv310"
echo "(Python 3.10)里装不了——这部分需要另外的conda环境，见下面第8步。"
echo ""
echo "=== 8. (可选)conda环境：配体力场参数化，真实跑蛋白-配体复合物MD需要 ==="
echo "本仓库复用一个已存在的conda环境(这里示例用名字'bio'，换成你自己的环境名也行)："
echo "  mamba install -n bio -c conda-forge ambertools openff-toolkit openmm openmmforcefields pdbfixer -y"
echo "  (用conda classic solver可能要卡在Solving environment几十分钟，强烈建议先装mamba："
echo "   conda install -n base -c conda-forge mamba -y，再用mamba install，几分钟内能解完)"
echo "  真实AmberTools的antechamber二进制装在该环境的bin目录下，调用前必须把这个目录"
echo "  加进PATH(PATH=\"/path/to/envs/bio/bin:\$PATH\")，否则openmmforcefields内部调用"
echo "  antechamber会报'command not found'。"
echo "  验证: PATH=\"/path/to/envs/bio/bin:\$PATH\" /path/to/envs/bio/bin/python -c \\"
echo "    \"from core.md_stability import run_protein_ligand_complex_md; ...\""
