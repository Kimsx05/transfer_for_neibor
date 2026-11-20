#!/bin/bash
set -eu
a=../cell_annotation_v2
export PYTHONPATH="../segmentation_restart_v1/packages:$a/.gpu_packages"
export LD_LIBRARY_PATH=$(printf '%s:' "$a"/.gpu_packages/nvidia/*/lib)
export PATH="$a/.gpu_packages/nvidia/cuda_nvcc/bin:$PATH"
exec ../cell_annotation_freeze_v1/.venv/bin/python code/detect.py "$@"
