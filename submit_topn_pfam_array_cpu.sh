#!/bin/bash
#$ -cwd
#$ -o logs/topn_pfam_cpu.$JOB_ID.$TASK_ID.out
#$ -j y
#$ -l h_data=8G,h_rt=12:00:00,highp
#$ -pe shared 8
#$ -t 1-3
#$ -M $USER@ucla.edu
#$ -m bea
#
# CPU SGE array job: run foldseek_topn_pfam.py once per Pfam-subset chunk.
#
# CPU twin of submit_topn_pfam_array.sh. Searches with --gpu 0 against the plain
# merged target DB, so it needs no padded GPU database. Use this when GPU nodes /
# the padded DB are unavailable. Each task N processes pfam_chunk_NN.tsv (from
# split_pfam_hits.py); outputs land under OUTPUT_BASE/chunk_NN/.
#
# Total memory = h_data x cores = 8G x 8 = 64G.
#
# Before submitting: build the chunks and set -t to match the chunk count, e.g.
#   uv run python split_pfam_hits.py "$SCRATCH/all_pfam_top50_hits.tsv" \
#       --num-chunks 3 --max-pfams 6 --out-dir "$SCRATCH/pfam_chunks_test"
#   qsub submit_topn_pfam_array_cpu.sh

# --- Environment (source before `set -u`; /etc/bashrc uses an unbound PS1) ---
. /u/local/Modules/default/init/modules.sh
source ~/.bashrc
ulimit -c 0

set -euo pipefail

# --- Configuration --------------------------------------------------------
TARGET_DB="$SCRATCH/all_3dis_fully_merged_2026-05-06/mergedDB"   # plain DB; no padding needed for CPU
CHUNKS_DIR="$SCRATCH/pfam_chunks_test"
OUTPUT_BASE="$SCRATCH/foldseek_topn_pfam_array_test_cpu"
FOLDSEEK_BIN="$HOME/prostT5-runner/foldseek/bin/foldseek"
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

echo "=== Task $SGE_TASK_ID (CPU) ==="
echo "=== Chunk: $CHUNK ==="
echo "=== Target DB: $TARGET_DB ==="
echo "=== Output: $TASK_OUTPUT ==="

# --- Run ------------------------------------------------------------------
uv run python -u ~/prostT5-runner/foldseek_topn_pfam.py \
    "$CHUNK" \
    "$TARGET_DB" \
    "$TASK_OUTPUT" \
    --foldseek "$FOLDSEEK_BIN" \
    --gpu 0 \
    --threads "$THREADS" \
    --split-memory-limit "$SPLIT_MEMORY_LIMIT" \
    --progress-interval 200

echo "=== Task $SGE_TASK_ID done ==="
