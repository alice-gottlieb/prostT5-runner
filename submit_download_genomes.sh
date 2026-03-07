#!/bin/bash
#$ -cwd
#$ -o logs/download_genomes.$JOB_ID.$TASK_ID.out
#$ -j y
#$ -l h_data=8G,h_rt=4:00:00
#$ -pe shared 2
#$ -t 1-5
#$ -M $USER@ucla.edu
#$ -m bea
#
# ---------------------------------------------------------------------------
# Job array for downloading bacterial reference genomes from NCBI using
# dehydrated downloads for speed. Each task downloads a chunk of accessions.
#
# Before submitting:
#   1. Fetch accessions and create chunks:
#        python fetch_bacterial_accessions.py --num-genomes 100 --chunks 5
#      Adjust --num-genomes (or use 'all') and --chunks as needed.
#      Update -t above to match the number of chunks.
#
#   2. Make sure NCBI datasets CLI is available (or install via conda/mamba):
#        conda install -c conda-forge ncbi-datasets-cli
#
#   3. Submit:
#        qsub submit_download_genomes.sh
#
#   4. After all tasks finish, genomes are in:
#        $SCRATCH/ncbi_genomes/task_NN/ncbi_dataset/data/GCF_*/
# ---------------------------------------------------------------------------

# ---- Configuration (edit these) ----
CHUNKS_DIR="$HOME/prostT5-runner/download_chunks"
OUTPUT_BASE="$SCRATCH/ncbi_genomes"
API_KEY_FILE="$HOME/prostT5-runner/ncbi_key.txt"
BATCH_SIZE=50  # accessions per datasets call (avoids URL length limits)
# ------------------------------------

# Source module system and user profile
. /u/local/Modules/default/init/modules.sh
source ~/.bashrc

TASK_ID_PAD=$(printf "%02d" $SGE_TASK_ID)
CHUNK_FILE="${CHUNKS_DIR}/chunk_${TASK_ID_PAD}.txt"

if [ ! -f "$CHUNK_FILE" ]; then
    echo "ERROR: Chunk file not found: $CHUNK_FILE"
    exit 1
fi

API_KEY=$(cat "$API_KEY_FILE")
TASK_OUTPUT="${OUTPUT_BASE}/task_${TASK_ID_PAD}"
mkdir -p "$TASK_OUTPUT" logs

NUM_ACCESSIONS=$(grep -cv '^\s*$\|^#' "$CHUNK_FILE")
echo "=== Task $SGE_TASK_ID: Downloading $NUM_ACCESSIONS genomes from $CHUNK_FILE ==="
echo "=== Started: $(date) ==="

# Read accessions into array
mapfile -t ALL_ACCESSIONS < <(grep -v '^\s*$\|^#' "$CHUNK_FILE")

# Process in batches
TOTAL=${#ALL_ACCESSIONS[@]}
BATCH_NUM=0

for ((i = 0; i < TOTAL; i += BATCH_SIZE)); do
    BATCH_NUM=$((BATCH_NUM + 1))
    BATCH=("${ALL_ACCESSIONS[@]:i:BATCH_SIZE}")
    ACC_LIST=$(IFS=,; echo "${BATCH[*]}")

    BATCH_DIR="${TASK_OUTPUT}/batch_${BATCH_NUM}"
    mkdir -p "$BATCH_DIR"

    echo "--- Batch $BATCH_NUM: ${#BATCH[@]} accessions ($(date)) ---"

    # Step 1: Dehydrated download (metadata + manifest only, very fast)
    datasets download genome accession $ACC_LIST \
        --api-key "$API_KEY" \
        --include protein \
        --dehydrated \
        --filename "${BATCH_DIR}/dehydrated.zip"

    if [ $? -ne 0 ]; then
        echo "ERROR: Dehydrated download failed for batch $BATCH_NUM"
        continue
    fi

    # Step 2: Unzip the dehydrated package
    unzip -o -q "${BATCH_DIR}/dehydrated.zip" -d "${BATCH_DIR}"
    rm -f "${BATCH_DIR}/dehydrated.zip"

    # Step 3: Rehydrate (fetch actual sequence files in parallel)
    datasets rehydrate --directory "${BATCH_DIR}" --api-key "$API_KEY"

    if [ $? -ne 0 ]; then
        echo "ERROR: Rehydration failed for batch $BATCH_NUM"
        continue
    fi

    echo "--- Batch $BATCH_NUM complete ---"
done

echo "=== Task $SGE_TASK_ID finished at $(date) ==="
echo "=== Output: $TASK_OUTPUT ==="
