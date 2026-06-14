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
PRECISIONS="${PRECISIONS:-fp16 int8}"   # which engines to build/validate

WORK="$SCRATCH/prostt5-trt"
REPO="$HOME/prostT5-runner"
export UV_CACHE_DIR="$WORK/uv-cache"
export HF_HOME="$WORK/hf-cache"

ONNX="$WORK/prostt5_3di.onnx"

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

if [ "$DO_EXPORT" = "1" ] || [ ! -f "$ONNX" ]; then
    echo "=== Exporting ONNX ==="
    python export_onnx.py \
        --cnn-ckpt "$WORK/cnn_head.pt" \
        --cache-dir "$HF_HOME" \
        --out "$ONNX"
fi

for PREC in $PRECISIONS; do
    ENGINE="$WORK/engines/prostt5_3di.${GPU_LABEL}.${PREC}.engine"
    echo "=== Building $PREC engine for $GPU_LABEL ==="
    BUILD_ARGS=(--onnx "$ONNX" --engine "$ENGINE"
                --min-len 16 --opt-len 384 --max-len "$MAX_LEN"
                --fp16 --workspace-gb "$WORKSPACE_GB")
    if [ "$PREC" = "int8" ]; then
        BUILD_ARGS+=(--int8
                     --calib-fasta "$WORK/calib_seqs.fasta"
                     --cache-dir "$HF_HOME"
                     --calib-cache "$WORK/engines/calib.${GPU_LABEL}.cache")
    fi
    python build_engine.py "${BUILD_ARGS[@]}"

    echo "=== Validating + benchmarking $PREC ($GPU_LABEL) ==="
    python validate_and_bench.py \
        --engine "$ENGINE" \
        --cnn-ckpt "$WORK/cnn_head.pt" \
        --cache-dir "$HF_HOME" \
        --fasta "$WORK/test_seqs.fasta" \
        --runs 5
done

echo "=== DONE for $GPU_LABEL ==="
