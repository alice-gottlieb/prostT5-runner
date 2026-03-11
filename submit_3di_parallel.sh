#!/bin/bash
#$ -cwd
#$ -o logs/all_3di/new_chunk_20_3/all_3di.$JOB_ID.$TASK_ID.out
#$ -j y
#$ -l gpu,H200,cuda=1,h_data=250G,h_rt=6:00:00,gpu_mem=24G
#$ -pe shared 4
#$ -t 450-850:4
#$ -M $USER@ucla.edu
#$ -m ea
#
# ---------------------------------------------------------------------------
# Single-job script that runs 8 batch_3di_foldseek.py processes in parallel
# (2 threads each, 16 cores total, CPU-only).
#
# Each process handles one chunk file. Edit NUM_CHUNKS and CHUNKS_DIR below
# to match your setup.
#
# Before submitting:
#   1. Split your accessions file into chunks:
#        python split_accessions.py reference_genomes.txt --chunks 8
#      This creates chunks/chunk_01.txt ... chunks/chunk_08.txt
#
#   2. Make sure foldseek and ProstT5 weights are available:
#        - foldseek binary: set FOLDSEEK_BIN below
#        - ProstT5 weights: set PROSTT5_WEIGHTS below (the .gguf file)
#
#   3. Submit:
#        qsub submit_3di_parallel.sh
#
#   4. After all chunks finish, combine 3Di FASTAs:
#        cat results_3di/chunk_*/foldseek_db/all_sequences_3di.fasta \
#            > results_3di/combined_3di.fasta
# ---------------------------------------------------------------------------

# ---- Configuration (edit these) ----
# FOLDSEEK_BIN="$HOME/prostT5-runner/foldseek/bin/foldseek"
# PROSTT5_WEIGHTS="$HOME/prostT5-runner/prostt5_weights/prostt5/prostt5-f16.gguf"
# NUM_CHUNKS=8
# THREADS_PER_JOB=2
# OUTPUT_BASE="$SCRATCH/compression_test"
# CHUNKS_DIR="$HOME/prostT5-runner/qsub_test"
FOLDSEEK_BIN="$HOME/prostT5-runner/foldseek/bin/foldseek"
PROSTT5_WEIGHTS="$HOME/prostT5-runner/prostt5_weights/prostt5/prostt5-f16.gguf"
THREADS_PER_JOB=2
NUM_CHUNKS=4
OUTPUT_BASE="$SCRATCH/all_3dis/new_chunks_20"
CHUNKS_DIR="$SCRATCH/new_chunks_20"
FASTA_DIR="$SCRATCH/ncbi_genomes"
# ------------------------------------

START_CHUNK=$SGE_TASK_ID
END_CHUNK=$((START_CHUNK + NUM_CHUNKS - 1))

echo "=== Task $SGE_TASK_ID: Processing chunks $START_CHUNK-$END_CHUNK in parallel ($THREADS_PER_JOB threads each) ==="
echo "=== Started: $(date) ==="

# Create logs dir
mkdir -p logs

# Disable core dumps
ulimit -c 0

# Source module system and user profile
. /u/local/Modules/default/init/modules.sh
source ~/.bashrc

# ---------------------------------------------------------------------------
# Launch one batch_3di_foldseek.py per chunk in parallel
# ---------------------------------------------------------------------------
PIDS=()
CHUNK_IDS=()

for ((i = START_CHUNK; i <= END_CHUNK; i++)); do
    CHUNK_PAD=$(printf "%02d" $i)
    CHUNK_FILE="${CHUNKS_DIR}/chunk_${CHUNK_PAD}.txt"
    CHUNK_OUTPUT="${OUTPUT_BASE}/chunk_${CHUNK_PAD}"

    if [ ! -f "$CHUNK_FILE" ]; then
        echo "  WARNING: Chunk file not found, skipping: $CHUNK_FILE"
        continue
    fi

    NUM_ACCESSIONS=$(grep -cv '^\s*$\|^#' "$CHUNK_FILE" || true)
    echo "  Launching chunk $CHUNK_PAD ($NUM_ACCESSIONS accessions) -> $CHUNK_OUTPUT"
    mkdir -p "$CHUNK_OUTPUT"

    uv run python -u ~/prostT5-runner/batch_3di_foldseek.py "$CHUNK_FILE" \
        --output-dir "$CHUNK_OUTPUT" \
        --foldseek-path "$FOLDSEEK_BIN" \
        --prostt5-weights "$PROSTT5_WEIGHTS" \
        --gpu \
        --threads $THREADS_PER_JOB \
        --fasta-dir "$FASTA_DIR" \
        --skip-foldseek \
        > "${CHUNK_OUTPUT}/run.log" 2>&1 &

    PIDS+=($!)
    CHUNK_IDS+=($CHUNK_PAD)
done

echo "Launched ${#PIDS[@]} parallel jobs. Waiting..."

# ---------------------------------------------------------------------------
# Wait for all jobs and collect exit codes
# ---------------------------------------------------------------------------
FAILED=0
for idx in "${!PIDS[@]}"; do
    wait "${PIDS[$idx]}"
    ec=$?
    if [ $ec -ne 0 ]; then
        echo "  FAILED: chunk ${CHUNK_IDS[$idx]} (exit code $ec)"
        ((FAILED++))
    else
        echo "  OK: chunk ${CHUNK_IDS[$idx]}"
    fi
done

echo "=== All jobs finished at $(date) ==="
echo "=== Failed: $FAILED / ${#PIDS[@]} ==="

# ---------------------------------------------------------------------------
# Combine 3Di FASTAs from all chunks
# ---------------------------------------------------------------------------
# COMBINED="${OUTPUT_BASE}/combined_3di.fasta"
# cat "${OUTPUT_BASE}"/chunk_*/foldseek_db/all_sequences_3di.fasta \
#     > "$COMBINED" 2>/dev/null

# if [ -f "$COMBINED" ]; then
#     TOTAL_3DI=$(grep -c '^>' "$COMBINED" || echo 0)
#     echo "=== Combined $TOTAL_3DI 3Di sequences into $COMBINED ==="
# fi

# # ---------------------------------------------------------------------------
# # Compress each chunk output
# # ---------------------------------------------------------------------------
# if [ $FAILED -eq 0 ]; then
#     for ((i = START_CHUNK; i <= END_CHUNK; i++)); do
#         CHUNK_PAD=$(printf "%04d" $i)
#         CHUNK_OUTPUT="${OUTPUT_BASE}/chunk_${CHUNK_PAD}"
#         if [ -d "$CHUNK_OUTPUT" ]; then
#             tar -czf "${CHUNK_OUTPUT}.tar.gz" -C "$(dirname "$CHUNK_OUTPUT")" "$(basename "$CHUNK_OUTPUT")"
#             # rm -rf "$CHUNK_OUTPUT"
#             echo "  Compressed chunk_${CHUNK_PAD}"
#         fi
#     done
# fi

# Exit with failure if any job failed
if [ $FAILED -gt 0 ]; then
    echo "ERROR: $FAILED jobs failed. Check individual run.log files."
    exit 1
fi

exit 0
