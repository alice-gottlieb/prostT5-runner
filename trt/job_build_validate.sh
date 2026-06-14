#!/bin/bash
# Build a ProstT5 3Di TensorRT engine on the current GPU node, then validate +
# benchmark it against the FP32 PyTorch reference.
#
# Pass the GPU resource and label on the qsub command line, e.g.:
#   qsub -l gpu,RTX2080Ti,cuda=1,h_rt=3:00:00,h_data=32G \
#        -v GPU_LABEL=RTX2080Ti,DO_EXPORT=1 trt/job_build_validate.sh
#
# DO_EXPORT=1 also exports the (architecture-independent) ONNX first; run it on
# the first job only.
#$ -cwd
#$ -N prostt5_trt
#$ -o /u/scratch/a/aliceg/logs/misc/$JOB_ID.out
#$ -j y
#$ -M $USER@ucla.edu
#$ -m ea
set -euo pipefail

GPU_LABEL="${GPU_LABEL:-unknown}"
DO_EXPORT="${DO_EXPORT:-0}"
MAX_LEN="${MAX_LEN:-2048}"          # lower for small-VRAM GPUs (RTX2080Ti=10GB)
WORKSPACE_GB="${WORKSPACE_GB:-8}"   # TRT build workspace; keep < VRAM

WORK="$SCRATCH/prostt5-trt"
REPO="$HOME/prostT5-runner"
export UV_CACHE_DIR="$WORK/uv-cache"
export HF_HOME="$WORK/hf-cache"

# FP16 ONNX (T5LayerNorm keeps FP32 variance). A strongly-typed engine built
# from this is numerically stable; a plain FP16 build of the FP32 ONNX is not
# (RMSNorm overflows in FP16 -> zeroed embeddings -> constant 3Di).
# Own subdir so its external-data weight files don't collide with the FP32 ONNX.
ONNX_FP16="$WORK/onnx_fp16/prostt5_3di.fp16.onnx"
ENGINE="$WORK/engines/prostt5_3di.${GPU_LABEL}.fp16.engine"
mkdir -p "$WORK/onnx_fp16"

# Keep an audit copy of this script alongside the outputs.
cp "$REPO/trt/job_build_validate.sh" "$WORK/job_build_validate.${GPU_LABEL}.sh"

# No `module load cuda` needed: the torch + tensorrt-cu12 wheels bundle their own
# CUDA runtime libs; the node driver provides libcuda. (`module` isn't even
# initialised in the batch shell.)
source "$WORK/venv/bin/activate"
cd "$REPO/trt"

echo "=== Node / GPU ==="
hostname
nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader
python -c "import torch, tensorrt as trt; print('torch', torch.__version__, '| TensorRT', trt.__version__, '| CUDA avail', torch.cuda.is_available())"

if [ "$DO_EXPORT" = "1" ] || [ ! -f "$ONNX_FP16" ]; then
    echo "=== Exporting FP16 ONNX (GPU) ==="
    python export_onnx.py --fp16 \
        --cnn-ckpt "$WORK/cnn_head.pt" \
        --cache-dir "$HF_HOME" \
        --out "$ONNX_FP16"
fi

echo "=== Building strongly-typed FP16 engine for $GPU_LABEL ==="
python build_engine.py \
    --onnx "$ONNX_FP16" --engine "$ENGINE" --strongly-typed \
    --min-len 16 --opt-len 384 --max-len "$MAX_LEN" \
    --workspace-gb "$WORKSPACE_GB"

echo "=== Validating + benchmarking FP16 ($GPU_LABEL) ==="
python validate_and_bench.py \
    --engine "$ENGINE" \
    --cnn-ckpt "$WORK/cnn_head.pt" \
    --cache-dir "$HF_HOME" \
    --fasta "$WORK/test_seqs.fasta" \
    --runs 5

echo "=== DONE for $GPU_LABEL ==="
