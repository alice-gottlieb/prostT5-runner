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

# Build an FP32 engine to isolate "ONNX graph broken" from "FP16 precision".
# (onnxruntime has no cp312/glibc-2.17 wheel on these el7 nodes, so we use the
# TRT stack itself as the oracle.)
FP32_ENGINE="$WORK/engines/diag_fp32.engine"
echo "=== Building FP32 diagnostic engine ==="
python build_engine.py \
    --onnx "$WORK/prostt5_3di.onnx" \
    --engine "$FP32_ENGINE" \
    --min-len 16 --opt-len 256 --max-len 512 --workspace-gb 8

echo "=== diag: PyTorch FP32 vs FP32-engine vs FP16-engine ==="
echo "--- against FP32 engine ---"
python diag_compare.py \
    --engine "$FP32_ENGINE" \
    --onnx "$WORK/prostt5_3di.onnx" \
    --cnn-ckpt "$WORK/cnn_head.pt" \
    --cache-dir "$HF_HOME"
echo "--- against FP16 engine ---"
python diag_compare.py \
    --engine "$WORK/engines/prostt5_3di.A100.fp16.engine" \
    --onnx "$WORK/prostt5_3di.onnx" \
    --cnn-ckpt "$WORK/cnn_head.pt" \
    --cache-dir "$HF_HOME"
