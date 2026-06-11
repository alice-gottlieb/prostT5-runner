#!/bin/bash
#$ -cwd
#$ -o logs/topn_pfam.$JOB_ID.$TASK_ID.out
#$ -j y
#$ -l gpu,RTX2080Ti,cuda=1,h_rt=12:00:00,h_data=16G
#$ -t 1-3
#$ -M $USER@ucla.edu
#$ -m bea
#
# GPU SGE array job: run foldseek_topn_pfam.py once per Pfam-subset chunk.
#
# Each array task N foldseek-searches the 3Di codes of its chunk
# (pfam_chunk_NN.tsv, produced by split_pfam_hits.py) against the merged target
# 3Di DB on the GPU. Outputs land under OUTPUT_BASE/chunk_NN/.
#
# Before submitting: build the chunks and set -t to match the chunk count, e.g.
#   uv run python split_pfam_hits.py "$SCRATCH/all_pfam_top50_hits.tsv" \
#       --num-chunks 3 --max-pfams 6 --out-dir "$SCRATCH/pfam_chunks_test"
#   qsub submit_topn_pfam_array.sh

# --- Environment ----------------------------------------------------------
# Source the module system and user profile BEFORE `set -u`: /etc/bashrc
# references an unbound PS1, which would abort the script under `set -u`.
. /u/local/Modules/default/init/modules.sh
source ~/.bashrc
ulimit -c 0

set -euo pipefail

# --- Configuration --------------------------------------------------------
TARGET_DB="$SCRATCH/all_3dis_fully_merged_2026-05-06/mergedDB"   # verify prefix before submit
CHUNKS_DIR="$SCRATCH/pfam_chunks_test"
OUTPUT_BASE="$SCRATCH/foldseek_topn_pfam_array_test"
FOLDSEEK_BIN="$HOME/foldseek/bin/foldseek"
THREADS=8
SPLIT_MEMORY_LIMIT="40G"

TASK_PAD=$(printf "%02d" "$SGE_TASK_ID")
CHUNK="$CHUNKS_DIR/pfam_chunk_${TASK_PAD}.tsv"
TASK_OUTPUT="$OUTPUT_BASE/chunk_${TASK_PAD}"
mkdir -p "$TASK_OUTPUT" logs

if [[ ! -f "$CHUNK" ]]; then
    echo "ERROR: chunk file not found: $CHUNK" >&2
    exit 1
fi

# Keep a copy of this submit script alongside the outputs for an audit trail.
cp "$(realpath "$0")" "$OUTPUT_BASE/$(basename "$0")"

# Pick the GPU with the most free memory.
export CUDA_VISIBLE_DEVICES=$(nvidia-smi --query-gpu=index,memory.free \
    --format=csv,noheader,nounits | sort -t',' -k2 -nr | head -1 | cut -d',' -f1)
echo "=== Task $SGE_TASK_ID using GPU: $CUDA_VISIBLE_DEVICES ==="
nvidia-smi --query-gpu=name,memory.total,memory.free,memory.used --format=csv

echo "=== Chunk: $CHUNK ==="
echo "=== Target DB: $TARGET_DB ==="
echo "=== Output: $TASK_OUTPUT ==="

# --- Run ------------------------------------------------------------------
uv run python -u ~/prostT5-runner/foldseek_topn_pfam.py \
    "$CHUNK" \
    "$TARGET_DB" \
    "$TASK_OUTPUT" \
    --foldseek "$FOLDSEEK_BIN" \
    --gpu 1 \
    --threads "$THREADS" \
    --split-memory-limit "$SPLIT_MEMORY_LIMIT" \
    --progress-interval 200

echo "=== Task $SGE_TASK_ID done ==="
