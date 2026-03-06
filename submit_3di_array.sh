#!/bin/bash
#$ -cwd
#$ -o logs/3di_array.$JOB_ID.$TASK_ID.out
#$ -j y
#$ -l gpu,A100,cuda=1,h_data=16G,h_rt=24:00:00
#$ -pe shared 8
#$ -t 1-12
#$ -M $USER@ucla.edu
#$ -m bea
#
# ---------------------------------------------------------------------------
# Job array for generating 3Di codes across 12 GPUs on Hoffman2.
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
#        qsub submit_3di_array.sh
#
#   4. After all tasks finish, combine 3Di FASTAs:
#        cat results_3di/task_*/foldseek_db/all_sequences_3di.fasta \
#            > results_3di/combined_3di.fasta
#
# To use different GPU types, change "A100" above to e.g. H100, L40S, V100.
# To request high-memory A100s: -l gpu,A100,gpu_mem=80G,cuda=1
# ---------------------------------------------------------------------------

# ---- Configuration (edit these) ----
FOLDSEEK_BIN="$HOME/foldseek/bin/foldseek"
PROSTT5_WEIGHTS="$HOME/prostT5-runner/prostt5_weights/prostt5/prostt5-f16.gguf"
THREADS=8
OUTPUT_BASE="$HOME/prostT5-runner/results_3di"
CHUNKS_DIR="$HOME/prostT5-runner/chunks"
# ------------------------------------

# Pad task ID to match chunk filenames (chunk_01.txt, chunk_02.txt, ...)
TASK_ID_PAD=$(printf "%02d" $SGE_TASK_ID)
CHUNK_FILE="${CHUNKS_DIR}/chunk_${TASK_ID_PAD}.txt"

if [ ! -f "$CHUNK_FILE" ]; then
    echo "ERROR: Chunk file not found: $CHUNK_FILE"
    exit 1
fi

NUM_ACCESSIONS=$(grep -cv '^\s*$\|^#' "$CHUNK_FILE")
echo "=== Task $SGE_TASK_ID: Processing $NUM_ACCESSIONS accessions from $CHUNK_FILE ==="
echo "=== GPU: $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader) ==="
echo "=== Started: $(date) ==="

# Create output directory for this task
TASK_OUTPUT="${OUTPUT_BASE}/task_${TASK_ID_PAD}"
mkdir -p "$TASK_OUTPUT" logs

# Load CUDA module
module load cuda/12.1

# Run batch_3di_foldseek.py with --skip-foldseek (3Di generation only)
python batch_3di_foldseek.py "$CHUNK_FILE" \
    --output-dir "$TASK_OUTPUT" \
    --foldseek-path "$FOLDSEEK_BIN" \
    --prostt5-weights "$PROSTT5_WEIGHTS" \
    --gpu \
    --threads $THREADS \
    --skip-foldseek

EXIT_CODE=$?

echo "=== Task $SGE_TASK_ID finished with exit code $EXIT_CODE at $(date) ==="
exit $EXIT_CODE
