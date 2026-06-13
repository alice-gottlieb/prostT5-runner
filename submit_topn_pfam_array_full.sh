#!/bin/bash
#$ -cwd
#$ -o logs/topn_pfam_full.$JOB_ID.$TASK_ID.out
#$ -j y
# H200 (141G VRAM), not RTX2080Ti (10G): the padded target DB's _ss file is 28G
# and must fit in GPU memory. On a 10G card foldseek's GPU batching either hangs
# or crashes (illegal memory access in libmarv/dbbatching.cuh); H200 holds the
# whole DB and has better queue turnaround than H100. Single core, generous host
# RAM: foldseek_topn_pfam.py also loads the 88M-pair target->genome map into a
# Python dict (~30G RSS) per task.
#$ -l gpu,H200,cuda=1,h_rt=8:00:00,h_data=48G
#$ -t 1-200
#$ -M $USER@ucla.edu
#$ -m bea
#
# GPU SGE array job: run foldseek_topn_pfam.py once per Pfam-subset chunk over
# the full top50 set (17,237 Pfams split into 200 chunks by split_pfam_hits.py).
# Each task searches its chunk's 3Di codes against the padded merged target DB
# on the GPU, then attributes hits to genomes via the full target->genome map.
#
# Pilot first:  qsub -t 1-4 submit_topn_pfam_array_full.sh   (measure runtime/mem)
# Then the rest: qsub -t 5-200 submit_topn_pfam_array_full.sh

. /u/local/Modules/default/init/modules.sh
source ~/.bashrc
ulimit -c 0
set -euo pipefail

# --- Configuration --------------------------------------------------------
# Padded GPU target DB (the _gpu_pad/ path from the old test script was purged;
# this _gpu_pad_24 build is the live one, rebuilt 2026-06-11).
TARGET_DB="$SCRATCH/all_3dis_fully_merged_2026-05-06_gpu_pad_24/mergedDB_pad"
CHUNKS_DIR="$SCRATCH/pfam_chunks"
OUTPUT_BASE="$SCRATCH/foldseek_topn_pfam_array_full"
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

# Keep a copy of this submit script alongside the outputs for an audit trail.
cp "$(realpath "$0")" "$OUTPUT_BASE/$(basename "$0")"

# Use the GPU that SGE allocated to this task. The scheduler sets
# CUDA_VISIBLE_DEVICES to the granted device (see `qstat -j`: resource map
# cuda=<node>=(N)). Do NOT override it from `nvidia-smi --query memory.free`:
# that ignores the allocation and makes concurrent array tasks on the same node
# pile onto one physical GPU, exhausting its VRAM (cudaMalloc / bad_alloc).
echo "=== Task $SGE_TASK_ID assigned CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-<UNSET>} ==="
nvidia-smi --query-gpu=index,name,memory.total,memory.free,memory.used --format=csv

echo "=== Chunk: $CHUNK ==="
echo "=== Target DB: $TARGET_DB ==="
echo "=== Map: $TARGET_GENOME_MAP ==="
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
