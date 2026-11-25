#!/usr/bin/env bash
set -euo pipefail
export TMPDIR="$PWD/tmp" PIP_CACHE_DIR="$PWD/tmp/pip" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MAKEFLAGS=-j1
py=/home/data/t070721/miniconda3/envs/compass_env/bin/python
"$py" -m venv cache/co310
p=cache/co310/bin/pip
site=cache/co310/lib/python3.10/site-packages
realpath ../../mp10_immune_endothelial_mechanisms_v1/venv/lib/python3.10/site-packages > "$site/project.pth"
export PIP_INDEX_URL=https://mirrors.ustc.edu.cn/pypi/simple
"$p" install --disable-pip-version-check 'Cython<3' wheel 'setuptools<71'
"$p" install --disable-pip-version-check 'numpy==1.26.4' 'pandas==1.5.3' 'scipy==1.11.4' 'matplotlib==3.6.3' 'scikit-learn==1.3.2' 'numba==0.59.1' 'anndata==0.9.2' 'scanpy==1.9.8' 'umap-learn==0.5.6' 'seaborn<0.14' 'igraph==0.11.8' 'louvain==0.8.2' goatools pysam 'ipython<9'
"$p" install --no-deps --no-build-isolation velocyto==0.17.17 ../mp10_srebf_nectin_virtual_perturbation_v1/cache/celloracle-0.18.0-py3-none-any.whl
export XDG_CONFIG_HOME="$PWD/tmp/xdgconfig" XDG_CACHE_HOME="$PWD/tmp/xdgcache"
"$p" install --no-build-isolation 'genomepy==0.16.4' 'gimmemotifs==0.17.2' 'ipywidgets==8.1.9'
cache/co310/bin/python -c 'import celloracle as co;print(co.__version__)'

