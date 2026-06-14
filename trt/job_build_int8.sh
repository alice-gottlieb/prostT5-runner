#!/bin/bash
# Build + validate an INT8 ProstT5 3Di engine on the current GPU.
#
# INT8 matmuls with a BF16 high-precision fallback (BF16 keeps the T5 RMSNorm
# stable, unlike FP16). BF16 is Ampere+; run this on A100 / H200, not the
# Turing RTX2080Ti. Built from the FP32 ONNX (the calibrator needs a
# builder-chosen-precision network, so it can't use the strongly-typed FP16 ONNX).
#
#   qsub -pe shared 4 -l gpu,A100,cuda=1,h_data=8G,h_rt=4:00:00 \
#        -v GPU_LABEL=A100,MAX_LEN=2048,WORKSPACE_GB=8 trt/job_build_int8.sh
#$ -cwd
#$ -N prostt5_int8
#$ -o /u/scratch/a/aliceg/logs/misc/$JOB_ID.out
#$ -j y
#$ -M $USER@ucla.edu
#$ -m ea
set -euo pipefail

GPU_LABEL="${GPU_LABEL:-unknown}"
MAX_LEN="${MAX_LEN:-2048}"
WORKSPACE_GB="${WORKSPACE_GB:-8}"

WORK="$SCRATCH/prostt5-trt"
REPO="$HOME/prostT5-runner"
export UV_CACHE_DIR="$WORK/uv-cache"
export HF_HOME="$WORK/hf-cache"

ONNX_FP32="$WORK/prostt5_3di.onnx"
ENGINE="$WORK/engines/prostt5_3di.${GPU_LABEL}.int8.engine"
cp "$REPO/trt/job_build_int8.sh" "$WORK/job_build_int8.${GPU_LABEL}.sh"

source "$WORK/venv/bin/activate"
cd "$REPO/trt"

echo "=== Node / GPU ==="; hostname
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

echo "=== Building INT8 (+BF16) engine for $GPU_LABEL ==="
python build_engine.py \
    --onnx "$ONNX_FP32" --engine "$ENGINE" --int8 --bf16 \
    --min-len 16 --opt-len 384 --max-len "$MAX_LEN" \
    --calib-fasta "$WORK/calib_seqs.fasta" \
    --cache-dir "$HF_HOME" \
    --calib-cache "$WORK/engines/calib.${GPU_LABEL}.cache" \
    --workspace-gb "$WORKSPACE_GB"

echo "=== Validating + benchmarking INT8 ($GPU_LABEL) ==="
python validate_and_bench.py \
    --engine "$ENGINE" \
    --cnn-ckpt "$WORK/cnn_head.pt" \
    --cache-dir "$HF_HOME" \
    --fasta "$WORK/test_seqs.fasta" \
    --runs 5

echo "=== DONE INT8 for $GPU_LABEL ==="
