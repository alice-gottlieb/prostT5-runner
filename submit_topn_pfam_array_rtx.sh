#!/bin/bash
#$ -cwd
#$ -o /u/scratch/a/aliceg/logs/pfam/topn_pfam_rtx.$JOB_ID.$TASK_ID.out
#$ -j y
# Production sharded-RTX array: run foldseek_topn_pfam.py per Pfam-subset chunk on
# RTX2080Ti (10G VRAM) by sharding ONLY the GPU prefilter. The full 29G padded DB
# does not fit a 10G card, but each ~4.6G padded shard does. Per task: GPU
# ungapped prefilter against each shard -> merge into a global top-1000 prefilter
# in the full-DB keyspace -> structurealign against the FULL DB (true full-DB
# E-values) -> aggregate. Reproduces the H200 full-DB run to ~98% of genome-count
# cells (strong hits and E-values exact; the ~2% residual is weak-tail GPU
# prefilter scoring differences). Peak VRAM ~5.5G, well under the card.
#
# RTX2080Ti is the most-available GPU (fast queue turnaround) and many array
# tasks run in parallel, unlike the single shared H200/H100 node.
#
# -pe shared 8 gives 8 real cores for the CPU structurealign; build_merged_prefilter
# plus the target->genome map need generous host RAM (h_data is PER core).
#$ -l gpu,RTX2080Ti,cuda=1,h_rt=12:00:00,h_data=9G
#$ -pe shared 8
#$ -t 1-200
#$ -M $USER@ucla.edu
#$ -m bea

. /u/local/Modules/default/init/modules.sh
source ~/.bashrc
ulimit -c 0
set -euo pipefail

# --- Configuration --------------------------------------------------------
SHARDS_DIR="$SCRATCH/rtx_shards_pad"               # padded shards for GPU prefilter
FULL_DB="$SCRATCH/all_3dis_fully_merged_2026-05-06/mergedDB"  # for align + E-values
CHUNKS_DIR="$SCRATCH/pfam_chunks"
OUTPUT_BASE="$SCRATCH/foldseek_topn_pfam_array_rtx"
TARGET_GENOME_MAP="$SCRATCH/target_genome_map_full/full_map.tsv"
FOLDSEEK_BIN="$HOME/prostT5-runner/foldseek/bin/foldseek"
THREADS=8

TASK_PAD=$(printf "%02d" "$SGE_TASK_ID")
CHUNK="$CHUNKS_DIR/pfam_chunk_${TASK_PAD}.tsv"
TASK_OUTPUT="$OUTPUT_BASE/chunk_${TASK_PAD}"
mkdir -p "$TASK_OUTPUT"

if [[ ! -f "$CHUNK" ]]; then
    echo "ERROR: chunk file not found: $CHUNK" >&2
    exit 1
fi

# Keep a copy of this submit script alongside the outputs for an audit trail.
cp "$(realpath "$0")" "$OUTPUT_BASE/$(basename "$0")"

# Use the GPU the scheduler granted; do NOT override from nvidia-smi --memory.free
# (that piles concurrent tasks onto one physical GPU and exhausts its VRAM).
echo "=== Task $SGE_TASK_ID CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-<UNSET>} ==="
nvidia-smi --query-gpu=index,name,memory.total,memory.free --format=csv

echo "=== Chunk: $CHUNK ==="
echo "=== Shards: $SHARDS_DIR ; Full DB: $FULL_DB ==="

# --- Run ------------------------------------------------------------------
uv run python -u ~/prostT5-runner/foldseek_topn_pfam.py \
    "$CHUNK" \
    "$FULL_DB" \
    "$TASK_OUTPUT" \
    --target-shards "$SHARDS_DIR" \
    --foldseek "$FOLDSEEK_BIN" \
    --threads "$THREADS" \
    --target-genome-map "$TARGET_GENOME_MAP" \
    --progress-interval 200

# Aggregate this chunk's hits into a Pfam x genome count matrix (chunks partition
# Pfams disjointly, so per-chunk matrices concatenate into the full matrix).
echo "=== Building Pfam x genome matrix for chunk $TASK_PAD ==="
uv run python -u ~/prostT5-runner/genome_pfam_matrix.py \
    "$TASK_OUTPUT/hits.tsv" \
    -o "$TASK_OUTPUT/genome_pfam_matrix.tsv"

echo "=== Task $SGE_TASK_ID done ==="
