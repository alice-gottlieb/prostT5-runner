#!/bin/bash
#$ -cwd
#$ -o logs/topn_pfam_rtx.$JOB_ID.$TASK_ID.out
#$ -j y
# RTX2080Ti pilot: searches against a SUBSET padded DB (~1/6 of targets, ~5G _ss)
# small enough to fit the ~11G RTX VRAM with no batching. Feasibility/perf test
# only -- hit counts are partial because the target DB is a subset.
#$ -l gpu,RTX2080Ti,cuda=1,h_rt=4:00:00,h_data=48G
#$ -t 1-4
#$ -M $USER@ucla.edu
#$ -m bea

. /u/local/Modules/default/init/modules.sh
source ~/.bashrc
ulimit -c 0
set -euo pipefail

# --- Configuration --------------------------------------------------------
# Subset padded DB from submit_build_rtx_subset_db.sh (fits RTX VRAM).
TARGET_DB="$SCRATCH/rtx_subset_pad/subDB_pad"
CHUNKS_DIR="$SCRATCH/pfam_chunks"
OUTPUT_BASE="$SCRATCH/foldseek_topn_pfam_array_rtx_pilot"
TARGET_GENOME_MAP="$SCRATCH/target_genome_map_full/full_map.tsv"
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

cp "$(realpath "$0")" "$OUTPUT_BASE/$(basename "$0")"

# Use the GPU that SGE allocated (it sets CUDA_VISIBLE_DEVICES). Do NOT override
# it from nvidia-smi -- that piles concurrent tasks onto one physical GPU.
echo "=== Task $SGE_TASK_ID assigned CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-<UNSET>} ==="
nvidia-smi --query-gpu=index,name,memory.total,memory.free,memory.used --format=csv

echo "=== Chunk: $CHUNK ==="
echo "=== Target DB (subset): $TARGET_DB ==="
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
    --target-genome-map "$TARGET_GENOME_MAP" \
    --progress-interval 200

# Aggregate this chunk's hits into a Pfam x genome count matrix. Chunks partition
# Pfams disjointly, so these per-chunk matrices concatenate (aligning genome
# columns) into the full matrix once every chunk is done.
echo "=== Building Pfam x genome matrix for chunk $TASK_PAD ==="
uv run python -u ~/prostT5-runner/genome_pfam_matrix.py \
    "$TASK_OUTPUT/hits.tsv" \
    -o "$TASK_OUTPUT/genome_pfam_matrix.tsv"

echo "=== Task $SGE_TASK_ID done ==="
