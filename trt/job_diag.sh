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

# onnxruntime (CPU) for the graph-isolation check.
uv pip install --python "$WORK/venv/bin/python" --only-binary :all: onnxruntime >/dev/null 2>&1 || \
    pip install onnxruntime >/dev/null 2>&1 || true

python diag_compare.py \
    --engine "$WORK/engines/prostt5_3di.A100.fp16.engine" \
    --onnx "$WORK/prostt5_3di.onnx" \
    --cnn-ckpt "$WORK/cnn_head.pt" \
    --cache-dir "$HF_HOME"
