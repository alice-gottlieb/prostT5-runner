#!/bin/bash
#$ -cwd
#$ -o /u/scratch/a/aliceg/test/$JOB_ID/logs/$JOB_ID.out
#$ -j y

##$ -l L40S=TRUE,cuda=1,h_data=1G,h_rt=6:00:00,h_vmem=INFINITY,require_gpu=TRUE
## 
##$ -l gpu,H200,cuda=1,h_data=32G,h_rt=5:00:00,gpu_mem=10G,h_vmem=32G
## ^ Success

#$ -l gpu,RTX2080Ti,cuda=1,h_data=32G,h_rt=0:30:00,gpu_mem=4G,h_vmem=16G

##$ -pe node 2
##$ -l gpu,gpu_model="RTX2080Ti|H200|A100|H100|L40S|A6000",cuda=1,h_data=32G,h_rt=6:00:00,gpu_mem=7G
##$ -l gpu,V100,cuda=1,h_data=16G,h_rt=6:00:00,gpu_mem=6G
##$ -pe shared 2
##$ -t 1-2

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
FOLDSEEK_BIN="$HOME/prostT5-runner/foldseek/bin/foldseek"
PROSTT5_WEIGHTS="$HOME/prostT5-runner/prostt5_weights/prostt5/prostt5-f16.gguf"
THREADS=1
# OUTPUT_BASE="$SCRATCH/all_3dis/4-11_final_chunk_20"
OUTPUT_BASE="$SCRATCH/test/rtxTest1/$JOB_ID"
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

# Pick the GPU with the most free memory
# export CUDA_VISIBLE_DEVICES=$(nvidia-smi --query-gpu=index,memory.free --format=csv,noheader,nounits | sort -t',' -k2 -nr | head -1 | cut -d',' -f1)
echo "=== CUDA Visible Devices: $CUDA_VISIBLE_DEVICES ==="

export CUDA_LAUNCH_BLOCKING=1

NUM_ACCESSIONS=$(grep -cv '^\s*$\|^#' "$CHUNK_FILE")
echo "=== Task $SGE_TASK_ID: Processing $NUM_ACCESSIONS accessions from $CHUNK_FILE ==="
echo "=== GPU: $(nvidia-smi --query-gpu=name,memory.total,memory.free,memory.used,driver_version,compute_cap,temperature.gpu,utilization.gpu,utilization.memory --format=csv) ==="
echo "=== Started: $(date +"%Y-%m-%d %H:%M:%S") ==="

# Create output directory for this task
# TASK_OUTPUT="${OUTPUT_BASE}/task_${TASK_ID_PAD}"
TASK_OUTPUT = "$OUTPUT_BASE/rtxTest1_output"
mkdir -p "$TASK_OUTPUT" logs

cp "$(realpath "$0")" "$OUTPUT_BASE/$(basename "$0")"

# Disable core dumps
ulimit -c 0

# Source module system and user profile
. /u/local/Modules/default/init/modules.sh
source ~/.bashrc

# Load CUDA module
# module load cuda/10.0

CHUNK_FILE = "$CHUNKS_DIR/rtxTest1.txt"

# Run batch_3di_foldseek.py with --skip-foldseek (3Di generation only)
uv run python -u ~/prostT5-runner/batch_3di_foldseek.py "$CHUNK_FILE" \
    --output-dir "$TASK_OUTPUT" \
    --foldseek-path "$FOLDSEEK_BIN" \
    --prostt5-weights "$PROSTT5_WEIGHTS" \
    --gpu \
    --threads $THREADS \
    --fasta-dir $FASTA_DIR \
    --skip-foldseek 
    # --ncbi-api-key $(cat ~/prostT5-runner/ncbi_key.txt)

EXIT_CODE=$?

echo "=== Task $SGE_TASK_ID finished with exit code $EXIT_CODE at $(date) ==="
exit $EXIT_CODE
