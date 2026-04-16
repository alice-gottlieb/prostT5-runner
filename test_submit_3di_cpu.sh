#!/bin/bash
#$ -cwd
#$ -o /u/scratch/a/aliceg/test/cpu/$JOB_ID/logs/$JOB_ID.out
#$ -j y
#$ -l h_data=8G,h_rt=23:00:00,highp
#$ -pe shared 8
##$ -t 294-295
#$ -M $USER@ucla.edu
#$ -m bea
#
# ---------------------------------------------------------------------------
# CPU job array for generating 3Di codes on Hoffman2.
#
# Before submitting:
#   1. Split your accessions file into 12 chunks:
#        python split_accessions.py reference_genomes.txt --chunks 12
#      This creates chunks/chunk_01.txt ... chunks/chunk_12.txt
#
#   2. Make sure foldseek and ProstT5 weights are available:
#        - foldseek binary: set FOLDSEEK_BIN below
#        - ProstT5 weights: set PROSTT5_WEIGHTS below (the .gguf file)
#
#   3. Submit:
#        qsub submit_3di_array_cpu.sh
#
#   4. After all tasks finish, combine 3Di FASTAs:
#        cat results_3di/task_*/foldseek_db/all_sequences_3di.fasta \
#            > results_3di/combined_3di.fasta
#
# Adjust -pe shared N and THREADS below to match the number of CPU cores you
# want per task. More cores per task = fewer tasks fit concurrently on the
# cluster, but each task finishes faster.
# ---------------------------------------------------------------------------

# ---- Configuration (edit these) ----
# FOLDSEEK_BIN="$HOME/prostT5-runner/foldseek/bin/foldseek"
# PROSTT5_WEIGHTS="$HOME/prostT5-runner/prostt5_weights/prostt5/prostt5-f16.gguf"
# THREADS=8
# OUTPUT_BASE="$SCRATCH/compression_test"
# CHUNKS_DIR="$HOME/prostT5-runner/qsub_test"
# FOLDSEEK_BIN="$HOME/prostT5-runner/foldseek/bin/foldseek"
# PROSTT5_WEIGHTS="$HOME/prostT5-runner/prostt5_weights/prostt5/prostt5-f16.gguf"
# THREADS=8
# OUTPUT_BASE="$SCRATCH/all_3dis/4-11_final_chunk_20"
# CHUNKS_DIR="$SCRATCH/final_chunk_20"
# FASTA_DIR="$SCRATCH/ncbi_genomes"

FOLDSEEK_BIN="$HOME/prostT5-runner/foldseek/bin/foldseek"
PROSTT5_WEIGHTS="$HOME/prostT5-runner/prostt5_weights/prostt5/prostt5-f16.gguf"
THREADS=8
# OUTPUT_BASE="$SCRATCH/all_3dis/4-11_final_chunk_20"
OUTPUT_BASE="$SCRATCH/test/cpuTest1/$JOB_ID"
# CHUNKS_DIR="$SCRATCH/final_chunk_20"
CHUNKS_DIR="$SCRATCH/test/rtxTest1_chunks"
FASTA_DIR="$SCRATCH/ncbi_genomes"
# ------------------------------------

# Pad task ID to match chunk filenames (chunk_01.txt, chunk_02.txt, ...)
# TASK_ID_PAD=$(printf "%02d" $SGE_TASK_ID)
# CHUNK_FILE="${CHUNKS_DIR}/chunk_${TASK_ID_PAD}.txt"

# if [ ! -f "$CHUNK_FILE" ]; then
#     echo "ERROR: Chunk file not found: $CHUNK_FILE"
#     exit 1
# fi

CHUNK_FILE = "$CHUNKS_DIR/rtxTest1.txt"

NUM_ACCESSIONS=$(grep -cv '^\s*$\|^#' "$CHUNK_FILE")
echo "=== Processing $NUM_ACCESSIONS accessions from $CHUNK_FILE ==="
echo "=== Host: $(hostname) ==="
echo "=== CPU: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | xargs) ==="
echo "=== Cores allocated: $THREADS ==="
echo "=== Memory allocated: $(qstat -j $JOB_ID 2>/dev/null | grep -E '^hard resource_list' | sed 's/.*h_data=\([^,]*\).*/\1/') per slot (total across $THREADS slots) ==="
echo "=== Memory (free -h): ==="
free -h
echo "=== Started: $(date) ==="

# Create output directory for this task
TASK_OUTPUT="${OUTPUT_BASE}/task_${TASK_ID_PAD}"
mkdir -p "$TASK_OUTPUT" logs

cp "$(realpath "$0")" "$OUTPUT_BASE/$(basename "$0").sh"

# Disable core dumps
ulimit -c 0

# Source module system and user profile
. /u/local/Modules/default/init/modules.sh
source ~/.bashrc

# Run batch_3di_foldseek.py with --skip-foldseek (3Di generation only, CPU)
uv run python -u ~/prostT5-runner/batch_3di_foldseek.py "$CHUNK_FILE" \
    --output-dir "$TASK_OUTPUT" \
    --foldseek-path "$FOLDSEEK_BIN" \
    --prostt5-weights "$PROSTT5_WEIGHTS" \
    --threads $THREADS \
    --skip-foldseek \
    --foldseek-args -v 3

EXIT_CODE=$?

echo "=== Task finished with exit code $EXIT_CODE at $(date) ==="
exit $EXIT_CODE

# # Compress output directory
# if [ $EXIT_CODE -eq 0 ]; then
#     tar -czf "${TASK_OUTPUT}.tar.gz" -C "$(dirname "$TASK_OUTPUT")" "$(basename "$TASK_OUTPUT")"
#     # rm -rf "$TASK_OUTPUT"
#     echo "=== Compressed output to ${TASK_OUTPUT}.tar.gz ==="
# fi

# exit $EXIT_CODE
