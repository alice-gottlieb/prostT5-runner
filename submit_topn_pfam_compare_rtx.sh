#!/bin/bash
#$ -cwd
#$ -o /u/scratch/a/aliceg/logs/misc/compare_rtx.$JOB_ID.out
#$ -j y
# Equivalence-test SHARDED run: the same chunk (pfam_chunk_05) the H200 reference
# used, searched against the 6 disjoint padded shards on an RTX2080Ti (10G VRAM)
# and merged back to full-DB semantics. Must reproduce the H200 output and stay
# within VRAM (the whole 29G DB does not fit; each ~5G shard does).
#$ -l gpu,RTX2080Ti,cuda=1,h_rt=12:00:00,h_data=9G
#$ -pe shared 8
#$ -M $USER@ucla.edu
#$ -m bea

. /u/local/Modules/default/init/modules.sh
source ~/.bashrc
ulimit -c 0
set -euo pipefail

# Which production chunk to run (qsub -v CHUNK_ID=01 to override; default 05).
CHUNK_ID="${CHUNK_ID:-05}"
SHARDS_DIR="$SCRATCH/rtx_shards_pad"
# Full (non-padded) merged DB: used for the structurealign + convertalis so the
# E-values are the true full-DB values. The padded shards are only for the GPU
# ungapped prefilter.
FULL_DB="$SCRATCH/all_3dis_fully_merged_2026-05-06/mergedDB"
CHUNK="$SCRATCH/pfam_chunks/pfam_chunk_${CHUNK_ID}.tsv"
TASK_OUTPUT="$SCRATCH/rtx_vs_h200_test/rtx_chunk${CHUNK_ID}"
TARGET_GENOME_MAP="$SCRATCH/target_genome_map_full/full_map.tsv"
FOLDSEEK_BIN="$HOME/prostT5-runner/foldseek/bin/foldseek"
mkdir -p "$TASK_OUTPUT"
cp "$(realpath "$0")" "$(dirname "$TASK_OUTPUT")/$(basename "$0")"

echo "=== Assigned CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-<UNSET>} ==="
nvidia-smi --query-gpu=index,name,memory.total,memory.free --format=csv

# Sample GPU memory every 5s so we can confirm the search stays within VRAM.
# Log the index too and filter to OUR allocated GPU -- a node has several RTX
# cards and we must not count another job's usage.
GPU_IDX="${CUDA_VISIBLE_DEVICES:-0}"
( while true; do
    nvidia-smi --query-gpu=index,memory.used --format=csv,noheader,nounits
    sleep 5
  done ) > "$TASK_OUTPUT/gpu_mem_trace.log" 2>&1 &
MON_PID=$!
trap 'kill "$MON_PID" 2>/dev/null || true' EXIT

uv run python -u ~/prostT5-runner/foldseek_topn_pfam.py \
    "$CHUNK" \
    "$FULL_DB" \
    "$TASK_OUTPUT" \
    --target-shards "$SHARDS_DIR" \
    --foldseek "$FOLDSEEK_BIN" \
    --threads 8 \
    --target-genome-map "$TARGET_GENOME_MAP" \
    --progress-interval 200

kill "$MON_PID" 2>/dev/null || true
echo "=== peak GPU memory.used (MiB) on our GPU index $GPU_IDX ==="
awk -F', *' -v g="$GPU_IDX" '$1==g {print $2}' "$TASK_OUTPUT/gpu_mem_trace.log" \
    | sort -n | tail -1
echo "=== RTX sharded run done ==="
