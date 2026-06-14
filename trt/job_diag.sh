#!/bin/bash
# Diagnostic GPU job: compare PyTorch vs onnxruntime vs TensorRT on one sequence.
# Submit on an A100 (where the engine already exists):
#   qsub -pe shared 4 -l gpu,A100,cuda=1,h_data=8G,h_rt=1:00:00 trt/job_diag.sh
#$ -cwd
#$ -N prostt5_diag
#$ -o /u/scratch/a/aliceg/logs/misc/$JOB_ID.out
#$ -j y
#$ -M $USER@ucla.edu
#$ -m ea
set -euo pipefail

WORK="$SCRATCH/prostt5-trt"
export UV_CACHE_DIR="$WORK/uv-cache"
export HF_HOME="$WORK/hf-cache"

source "$WORK/venv/bin/activate"
cd "$HOME/prostT5-runner/trt"

hostname
nvidia-smi --query-gpu=name --format=csv,noheader

# onnxruntime (CPU) for the graph-isolation check (install loudly so we see why
# if it fails). uv isn't on the batch PATH by default.
export PATH="$HOME/.local/bin:$PATH"
# onnxruntime <=1.18 still ships manylinux_2_17 wheels (the el7 nodes are glibc
# 2.17); newer versions need glibc 2.28.
uv pip install --python "$WORK/venv/bin/python" --only-binary :all: "onnxruntime==1.18.1"
python -c "import onnxruntime; print('onnxruntime', onnxruntime.__version__)"

python diag_compare.py \
    --engine "$WORK/engines/prostt5_3di.A100.fp16.engine" \
    --onnx "$WORK/prostt5_3di.onnx" \
    --cnn-ckpt "$WORK/cnn_head.pt" \
    --cache-dir "$HF_HOME"
