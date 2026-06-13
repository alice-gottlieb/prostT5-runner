#!/bin/bash
#$ -cwd
#$ -o /u/scratch/a/aliceg/logs/misc/probe_rtx_native_split.$JOB_ID.out
#$ -j y
# Phase A feasibility probe: can `foldseek search --gpu 1 --split 6 --split-mode 0`
# run the FULL 29G padded target DB on a 10G RTX2080Ti by splitting the target
# into ~5G chunks, instead of crashing in the GPU auto-batcher (dbbatching.cuh)?
# Tiny query set so it finishes fast; we only care whether it completes and how
# much VRAM each split uses.
#$ -l gpu,RTX2080Ti,cuda=1,h_rt=1:00:00,h_data=16G
#$ -M $USER@ucla.edu
#$ -m bea

. /u/local/Modules/default/init/modules.sh
source ~/.bashrc
ulimit -c 0
set -euo pipefail

TARGET_DB="$SCRATCH/all_3dis_fully_merged_2026-05-06_gpu_pad_24/mergedDB_pad"
FOLDSEEK_BIN="$HOME/prostT5-runner/foldseek/bin/foldseek"
OUTDIR="$SCRATCH/probe_rtx_native_split"
mkdir -p "$OUTDIR"

# Tiny query set: header + first 100 data rows of a real chunk.
TINY_INPUT="$OUTDIR/tiny_query.tsv"
head -n 101 "$SCRATCH/pfam_chunks/pfam_chunk_05.tsv" > "$TINY_INPUT"
echo "=== tiny query rows: $(($(wc -l < "$TINY_INPUT") - 1)) ==="

echo "=== Assigned CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-<UNSET>} ==="
nvidia-smi --query-gpu=index,name,memory.total,memory.free --format=csv

# Sample GPU memory every 5s in the background so we can confirm we stay <10G.
( while true; do
    date +%T
    nvidia-smi --query-gpu=memory.used,memory.free --format=csv,noheader
    sleep 5
  done ) > "$OUTDIR/gpu_mem_trace.log" 2>&1 &
MON_PID=$!

set +e
uv run python -u ~/prostT5-runner/foldseek_topn_pfam.py \
    "$TINY_INPUT" \
    "$TARGET_DB" \
    "$OUTDIR/out" \
    --foldseek "$FOLDSEEK_BIN" \
    --gpu 1 \
    --threads 8 \
    --split 6 \
    --split-mode 0 \
    --split-memory-limit 8G \
    --progress-interval 50
RC=$?
set -e
kill "$MON_PID" 2>/dev/null || true

echo "=== peak GPU memory.used during run ==="
sort -t',' -k1 -n "$OUTDIR/gpu_mem_trace.log" 2>/dev/null | grep MiB | tail -1 || true
echo "=== foldseek_topn_pfam exit code: $RC ==="
exit $RC
