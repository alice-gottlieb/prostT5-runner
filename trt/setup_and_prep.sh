#!/bin/bash
# One-time setup + data prep for the ProstT5 TensorRT work.
# Runs on the Hoffman LOGIN node (only does I/O: venv creation, pip, downloads).
# All artifacts live in $SCRATCH/prostt5-trt -- nothing is written to $HOME or
# the distil-prostt5 directory.
set -euo pipefail

WORK="$SCRATCH/prostt5-trt"
REPO="$HOME/prostT5-runner"
export UV_CACHE_DIR="$WORK/uv-cache"          # keep uv cache off $HOME
export HF_HOME="$WORK/hf-cache"               # our own HF cache (not distil's)

mkdir -p "$WORK" "$WORK/engines" "$WORK/uv-cache" "$WORK/hf-cache"
mkdir -p "$SCRATCH/logs/misc"

echo "=== Creating venv at $WORK/venv ==="
if [ ! -x "$WORK/venv/bin/python" ]; then
    uv venv --python 3.12 "$WORK/venv"
fi
source "$WORK/venv/bin/activate"
cd "$REPO/trt"   # so the python heredocs below can `import prostt5_3di`

echo "=== Installing deps ==="
uv pip install --python "$WORK/venv/bin/python" \
    torch==2.6.0 --index-url https://download.pytorch.org/whl/cu124
# TensorRT 10.x for CUDA 12 (the meta "tensorrt" pulls cu13/TRT11). The libs +
# bindings live on NVIDIA's index, so add it as an extra index.
# This CentOS-7 login node has GCC 4.8.5, too old to compile modern packages, so
# force prebuilt wheels for everything (uv then picks versions that ship wheels).
uv pip install --python "$WORK/venv/bin/python" --only-binary :all: \
    "transformers==5.0.0" "sentencepiece==0.2.0" protobuf "numpy==2.2.6" onnx

# TensorRT 10.x (cu12). The top-level shim is a tiny sdist that pulls the real
# libs/bindings wheels from NVIDIA's index, so it is installed separately
# (its build step succeeds even with the old toolchain).
uv pip install --python "$WORK/venv/bin/python" \
    --extra-index-url https://pypi.nvidia.com \
    "tensorrt-cu12>=10.6,<11"

echo "=== Versions ==="
# Don't import tensorrt here: the login node's memory ulimit can't mmap its large
# .so (works fine on the GPU compute nodes). Just confirm torch/transformers and
# that the tensorrt package is installed.
python -c "import torch, transformers; \
print('torch', torch.__version__); \
print('transformers', transformers.__version__)"
python -c "import importlib.metadata as m; print('tensorrt-cu12', m.version('tensorrt-cu12'))" || true

echo "=== Env ready. Work dir: $WORK ==="
echo "Next: submit trt/job_bootstrap.sh (qsub) to download weights + build the"
echo "test/calibration FASTAs, then trt/job_build_validate.sh per GPU arch."
ls -la "$WORK"
