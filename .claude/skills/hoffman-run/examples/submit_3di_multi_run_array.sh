#!/bin/bash
# ILLUSTRATIVE EXAMPLE — a real multi-run GPU array job. The resources, paths,
# and accessions below are aliceg-specific; copy this and adapt them before
# submitting. Some header comments (A100/L40S/V100, split_accessions.py) and the
# commented-out blocks are leftovers and don't match the actual H200 -l line.
# See array-jobs.md for the array mechanics and resource-spec.md for GPU flags.
#$ -cwd
#$ -o /u/scratch/a/aliceg/logs/all_3di/final9k_chunk20/final.multirun.h200.all_3di.$JOB_ID.$TASK_ID.out
#$ -j y
#$ -l gpu,H200,cuda=1,h_data=16G,h_rt=6:00:00,gpu_mem=16G,h_vmem=32G
#$ -pe shared 8
#$ -t 200-480:2
#$ -M $USER@ucla.edu
#$ -m ea
#
# ---------------------------------------------------------------------------
# Job array for generating 3Di codes across GPUs on Hoffman2.
#
# Each SGE task launches PARALLEL_RUNS foldseek runs concurrently on the same
# GPU, processing chunk_${SGE_TASK_ID}, chunk_${SGE_TASK_ID+1}, ...,
# chunk_${SGE_TASK_ID+PARALLEL_RUNS-1}. Set -t step size equal to PARALLEL_RUNS
# so that chunk ranges don't overlap between tasks (here: step 4 matches
# PARALLEL_RUNS=4).
#
# Before submitting:
#   1. Split your accessions file into chunks:
#        python split_accessions.py reference_genomes.txt --chunks N
#      This creates chunks/chunk_01.txt ... chunks/chunk_NN.txt
#
#   2. Make sure foldseek and ProstT5 weights are available:
#        - foldseek binary: set FOLDSEEK_BIN below
#        - ProstT5 weights: set PROSTT5_WEIGHTS below (the .gguf file)
#
#   3. Submit:
#        qsub submit_3di_multi_run_array.sh
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
THREADS=4
PARALLEL_RUNS=2
# OUTPUT_BASE="$SCRATCH/compression_test"
# CHUNKS_DIR="$HOME/prostT5-runner/qsub_test"
OUTPUT_BASE="$SCRATCH/all_3dis/final9k_chunk20"
CHUNKS_DIR="$SCRATCH/final9k_chunk20"
FASTA_DIR="$SCRATCH/ncbi_genomes"
LOG_BASE="/u/scratch/a/aliceg/logs/all_3di/final9k_chunk20"
# ------------------------------------

# Source module system and user profile
. /u/local/Modules/default/init/modules.sh
source ~/.bashrc

# Load CUDA module
module load cuda/12.3

# Disable core dumps
ulimit -c 0

# Pick the GPU with the most free memory (shared across all parallel runs)
export CUDA_VISIBLE_DEVICES=$(nvidia-smi --query-gpu=index,memory.free --format=csv,noheader,nounits | sort -t',' -k2 -nr | head -1 | cut -d',' -f1)
# export CUDA_VISIBLE_DEVICES=$(nvidia-smi --query-gpu=index,memory.free --format=csv,noheader,nounits | sort -t',' -k2 -nr | head -2 | cut -d',' -f1 | paste -sd',')
echo "=== Using GPU: $CUDA_VISIBLE_DEVICES ==="
echo "=== GPU: $(nvidia-smi --query-gpu=name,memory.total,memory.free,memory.used,driver_version,compute_cap,temperature.gpu,utilization.gpu,utilization.memory --format=csv) ==="
echo "=== CUDA version: $(nvidia-smi --query-gpu=driver_version --format=csv,noheader 2>/dev/null | head -1) (driver); $(nvidia-smi 2>/dev/null | grep -oP 'CUDA Version: \K[0-9.]+') (runtime) ==="
echo "=== SGE task $SGE_TASK_ID launching $PARALLEL_RUNS parallel runs: chunks $SGE_TASK_ID..$((SGE_TASK_ID + PARALLEL_RUNS - 1)) ==="
echo "=== Started: $(date) ==="

mkdir -p "$OUTPUT_BASE" "$LOG_BASE"
cp "$(realpath "$0")" "$OUTPUT_BASE/$(basename "$0")"

# # Disable core dumps
# ulimit -c 0

# # Source module system and user profile
# . /u/local/Modules/default/init/modules.sh
# source ~/.bashrc

# # Load CUDA module
# module load cuda/12.3

# ---- Launch PARALLEL_RUNS foldseek jobs in the background ----
PIDS=()
CHUNK_IDS=()
TASK_OUTPUTS=()

for (( i=0; i<PARALLEL_RUNS; i++ )); do
    CHUNK_ID=$((SGE_TASK_ID + i))
    CHUNK_ID_PAD=$(printf "%02d" "$CHUNK_ID")
    CHUNK_FILE="${CHUNKS_DIR}/chunk_${CHUNK_ID_PAD}.txt"

    if [ ! -f "$CHUNK_FILE" ]; then
        echo "WARN: Chunk file not found, skipping: $CHUNK_FILE"
        continue
    fi

    TASK_OUTPUT="${OUTPUT_BASE}/task_${CHUNK_ID_PAD}"
    RUN_LOG="${LOG_BASE}/final.multirun.h200.all_3di_${JOB_ID}.${SGE_TASK_ID}.chunk_${CHUNK_ID_PAD}.log"
    mkdir -p "$TASK_OUTPUT"

    NUM_ACCESSIONS=$(grep -cv '^\s*$\|^#' "$CHUNK_FILE")
    echo "=== [chunk $CHUNK_ID_PAD] $NUM_ACCESSIONS accessions from $CHUNK_FILE -> $TASK_OUTPUT (log: $RUN_LOG) ==="

    (
        uv run python -u ~/prostT5-runner/batch_3di_foldseek.py "$CHUNK_FILE" \
            --output-dir "$TASK_OUTPUT" \
            --foldseek-path "$FOLDSEEK_BIN" \
            --prostt5-weights "$PROSTT5_WEIGHTS" \
            --gpu \
            --threads $THREADS \
            --fasta-dir $FASTA_DIR \
            --skip-foldseek
    ) > "$RUN_LOG" 2>&1 &

    PIDS+=($!)
    CHUNK_IDS+=("$CHUNK_ID_PAD")
    TASK_OUTPUTS+=("$TASK_OUTPUT")
done

# ---- Progress bar helper (renders to main SGE log) ----
TOTAL_RUNS=${#PIDS[@]}
BAR_WIDTH=40
START_TS=$(date +%s)

render_progress() {
    local done=$1 total=$2
    local now=$(date +%s)
    local elapsed=$(( now - START_TS ))
    local filled=0 empty=$BAR_WIDTH
    if [ "$total" -gt 0 ]; then
        filled=$(( done * BAR_WIDTH / total ))
        empty=$(( BAR_WIDTH - filled ))
    fi
    local bar=""
    [ $filled -gt 0 ] && bar=$(printf '%*s' "$filled" '' | tr ' ' '#')
    [ $empty  -gt 0 ] && bar+=$(printf '%*s' "$empty"  '' | tr ' ' '-')
    local pct=0
    [ "$total" -gt 0 ] && pct=$(( done * 100 / total ))
    printf '[PROGRESS] [%s] %d/%d runs (%d%%) | elapsed %ds\n' \
        "$bar" "$done" "$total" "$pct" "$elapsed"
    printf '[PROGRESS] GPUs: %s | %s\n' \
        "$CUDA_VISIBLE_DEVICES" \
        "$(nvidia-smi --query-gpu=name,memory.free,memory.used,utilization.gpu,temperature.gpu --format=csv,noheader --id="$CUDA_VISIBLE_DEVICES" 2>/dev/null | tr '\n' ' | ')"
    printf '[PROGRESS] CUDA: %s (driver); %s (runtime)\n' \
        "$(nvidia-smi --query-gpu=driver_version --format=csv,noheader --id="$CUDA_VISIBLE_DEVICES" 2>/dev/null | head -1)" \
        "$(nvidia-smi 2>/dev/null | grep -oP 'CUDA Version: \K[0-9.]+')"
}

render_progress 0 "$TOTAL_RUNS"

# ---- Wait for all parallel runs to finish (with progress bar) ----
OVERALL_EXIT=0
COMPLETED=0
REMAINING=("${!PIDS[@]}")
POLL_INTERVAL=300

while [ ${#REMAINING[@]} -gt 0 ]; do
    sleep $POLL_INTERVAL
    STILL_RUNNING=()
    for idx in "${REMAINING[@]}"; do
        pid="${PIDS[$idx]}"
        chunk_id="${CHUNK_IDS[$idx]}"
        task_output="${TASK_OUTPUTS[$idx]}"

        if kill -0 "$pid" 2>/dev/null; then
            STILL_RUNNING+=("$idx")
            continue
        fi

        if wait "$pid" 2>/dev/null; then
            RC=0
        else
            RC=$?
            OVERALL_EXIT=$RC
        fi
        COMPLETED=$(( COMPLETED + 1 ))
        echo "=== [chunk $chunk_id] finished with exit code $RC at $(date) ==="

        # if [ $RC -eq 0 ]; then
        #     tar -czf "${task_output}.tar.gz" -C "$(dirname "$task_output")" "$(basename "$task_output")"
        #     # rm -rf "$task_output"
        #     echo "=== [chunk $chunk_id] compressed output to ${task_output}.tar.gz ==="
        # fi

        render_progress "$COMPLETED" "$TOTAL_RUNS"
    done
    REMAINING=("${STILL_RUNNING[@]}")

    # Heartbeat while runs are still in flight
    if [ ${#REMAINING[@]} -gt 0 ]; then
        render_progress "$COMPLETED" "$TOTAL_RUNS"
    fi
done

echo "=== SGE task $SGE_TASK_ID finished all parallel runs with overall exit code $OVERALL_EXIT at $(date) ==="
exit $OVERALL_EXIT
