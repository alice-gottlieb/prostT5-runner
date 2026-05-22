#!/bin/bash
#$ -cwd
#$ -o /u/scratch/a/aliceg/logs/pfam/hmmer2-1_by_1.$JOB_ID.out
#$ -j y
#$ -l h_data=2G,h_rt=23:00:00,highp
#$ -pe shared 2
#$ -t 1-1000:10
#$ -M $USER@ucla.edu
#$ -m ea
. /u/local/Modules/default/init/modules.sh
source ~/.bashrc

GENOME_LIST="/u/scratch/a/aliceg/pfam_hits2-1by1/genomes_remaining.txt"
OUTPUT_DIR="/u/scratch/a/aliceg/pfam_hits2-1by1"
PFAM_HMM="/u/scratch/a/aliceg/pfam_data/Pfam-A.hmm"
CPUS=2
FILES_PER_TASK=10
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HMMSEARCH_SINGLE="hmmsearch_single.sh"

if [[ -z "${SGE_TASK_ID:-}" ]]; then
    echo "ERROR: SGE_TASK_ID is not set. This script should be run as an SGE array job." >&2
    exit 1
fi

if [[ ! "$SGE_TASK_ID" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: SGE_TASK_ID must be a positive integer: $SGE_TASK_ID" >&2
    exit 1
fi

if [[ ! -f "$GENOME_LIST" ]]; then
    echo "ERROR: genome list not found: $GENOME_LIST" >&2
    exit 1
fi

if [[ ! -x "$HMMSEARCH_SINGLE" ]]; then
    echo "ERROR: hmmsearch_single.sh not found or not executable: $HMMSEARCH_SINGLE" >&2
    exit 1
fi

START_LINE=$SGE_TASK_ID
END_LINE=$((SGE_TASK_ID + FILES_PER_TASK - 1))

echo "=== SGE task $SGE_TASK_ID processing genome list lines $START_LINE-$END_LINE ==="
echo "=== Genome list: $GENOME_LIST ==="
echo "=== Output dir: $OUTPUT_DIR ==="
echo "=== Started: $(date) ==="

mkdir -p "$OUTPUT_DIR"

for LINE_NUM in $(seq "$START_LINE" "$END_LINE"); do
    FAA_FILE=$(sed -n "${LINE_NUM}p" "$GENOME_LIST")

    if [[ -z "$FAA_FILE" ]]; then
        echo "WARN: line $LINE_NUM is empty or beyond end of $GENOME_LIST; skipping"
        continue
    fi

    echo "=== [line $LINE_NUM] hmmsearch_single.sh $FAA_FILE ==="
    if ! "$HMMSEARCH_SINGLE" "$FAA_FILE" "$OUTPUT_DIR" "$PFAM_HMM" "$CPUS"; then
        echo "WARN: hmmsearch failed for line $LINE_NUM: $FAA_FILE; continuing"
        continue
    fi
done

echo "=== SGE task $SGE_TASK_ID finished lines $START_LINE-$END_LINE at $(date) ==="
exit 0
