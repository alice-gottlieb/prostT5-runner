#!/bin/bash
#$ -cwd
#$ -o /u/scratch/a/aliceg/logs/misc/compare_h200.$JOB_ID.out
#$ -j y
# Equivalence-test REFERENCE run: one full production chunk (pfam_chunk_05)
# searched against the WHOLE 29G padded target DB on an H200 (141G VRAM, fits
# whole). This is the ground truth the sharded RTX run must reproduce.
#$ -l gpu,H200,cuda=1,h_rt=8:00:00,h_data=48G
#$ -M $USER@ucla.edu
#$ -m bea

. /u/local/Modules/default/init/modules.sh
source ~/.bashrc
ulimit -c 0
set -euo pipefail

TARGET_DB="$SCRATCH/all_3dis_fully_merged_2026-05-06_gpu_pad_24/mergedDB_pad"
CHUNK="$SCRATCH/pfam_chunks/pfam_chunk_05.tsv"
TASK_OUTPUT="$SCRATCH/rtx_vs_h200_test/h200_chunk05"
TARGET_GENOME_MAP="$SCRATCH/target_genome_map_full/full_map.tsv"
FOLDSEEK_BIN="$HOME/prostT5-runner/foldseek/bin/foldseek"
mkdir -p "$TASK_OUTPUT"
cp "$(realpath "$0")" "$(dirname "$TASK_OUTPUT")/$(basename "$0")"

echo "=== Assigned CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-<UNSET>} ==="
nvidia-smi --query-gpu=index,name,memory.total,memory.free --format=csv

uv run python -u ~/prostT5-runner/foldseek_topn_pfam.py \
    "$CHUNK" \
    "$TARGET_DB" \
    "$TASK_OUTPUT" \
    --foldseek "$FOLDSEEK_BIN" \
    --gpu 1 \
    --threads 8 \
    --split-memory-limit 40G \
    --target-genome-map "$TARGET_GENOME_MAP" \
    --progress-interval 200

echo "=== H200 reference done ==="
