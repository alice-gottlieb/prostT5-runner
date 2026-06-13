#!/bin/bash
#$ -cwd
#$ -o /u/scratch/a/aliceg/logs/misc/compare_rtx.$JOB_ID.out
#$ -j y
# Equivalence-test SHARDED run: the same chunk (pfam_chunk_05) the H200 reference
# used, searched against the 6 disjoint padded shards on an RTX2080Ti (10G VRAM)
# and merged back to full-DB semantics. Must reproduce the H200 output and stay
# within VRAM (the whole 29G DB does not fit; each ~5G shard does).
#$ -l gpu,RTX2080Ti,cuda=1,h_rt=8:00:00,h_data=48G
#$ -M $USER@ucla.edu
#$ -m bea

. /u/local/Modules/default/init/modules.sh
source ~/.bashrc
ulimit -c 0
set -euo pipefail

SHARDS_DIR="$SCRATCH/rtx_shards_pad"
CHUNK="$SCRATCH/pfam_chunks/pfam_chunk_05.tsv"
TASK_OUTPUT="$SCRATCH/rtx_vs_h200_test/rtx_chunk05"
TARGET_GENOME_MAP="$SCRATCH/target_genome_map_full/full_map.tsv"
FOLDSEEK_BIN="$HOME/prostT5-runner/foldseek/bin/foldseek"
mkdir -p "$TASK_OUTPUT"
cp "$(realpath "$0")" "$(dirname "$TASK_OUTPUT")/$(basename "$0")"

echo "=== Assigned CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-<UNSET>} ==="
nvidia-smi --query-gpu=index,name,memory.total,memory.free --format=csv

# Sample GPU memory every 5s so we can confirm the search stays within VRAM.
( while true; do
    nvidia-smi --query-gpu=memory.used --format=csv,noheader
    sleep 5
  done ) > "$TASK_OUTPUT/gpu_mem_trace.log" 2>&1 &
MON_PID=$!
trap 'kill "$MON_PID" 2>/dev/null || true' EXIT

uv run python -u ~/prostT5-runner/foldseek_topn_pfam.py \
    "$CHUNK" \
    "$TASK_OUTPUT" \
    --target-shards "$SHARDS_DIR" \
    --foldseek "$FOLDSEEK_BIN" \
    --gpu 1 \
    --threads 8 \
    --split-memory-limit 8G \
    --target-genome-map "$TARGET_GENOME_MAP" \
    --progress-interval 200

kill "$MON_PID" 2>/dev/null || true
echo "=== peak GPU memory.used (MiB) ==="
grep MiB "$TASK_OUTPUT/gpu_mem_trace.log" | sort -n | tail -1
echo "=== RTX sharded run done ==="
