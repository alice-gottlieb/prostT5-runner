#!/bin/bash
# Smoke-test the standalone predict_3di CLI end-to-end: AA FASTA -> 3Di FASTA,
# using a built engine. Pass GPU_LABEL to pick which engine.
#   qsub -pe shared 4 -l gpu,A100,cuda=1,h_data=8G,h_rt=1:00:00 \
#        -v GPU_LABEL=A100 trt/job_predict_smoke.sh
#$ -cwd
#$ -N prostt5_pred
#$ -o /u/scratch/a/aliceg/logs/misc/$JOB_ID.out
#$ -j y
#$ -M $USER@ucla.edu
#$ -m ea
set -euo pipefail

GPU_LABEL="${GPU_LABEL:-A100}"
PREC="${PREC:-fp16}"
WORK="$SCRATCH/prostt5-trt"
export HF_HOME="$WORK/hf-cache"

source "$WORK/venv/bin/activate"
cd "$HOME/prostT5-runner/trt"

ENGINE="$WORK/engines/prostt5_3di.${GPU_LABEL}.${PREC}.engine"
OUT="$WORK/smoke_3di.${GPU_LABEL}.${PREC}.fasta"

hostname
nvidia-smi --query-gpu=name --format=csv,noheader

python predict_3di.py \
    --engine "$ENGINE" \
    --cache-dir "$HF_HOME" \
    --in "$WORK/test_seqs.fasta" \
    --out "$OUT"

echo "=== first 6 lines of the 3Di FASTA ==="
head -n 6 "$OUT"
echo "=== sanity: distinct 3Di letters used (should be many, not just D) ==="
grep -v ">" "$OUT" | head -5 | fold -w1 | sort -u | tr -d "\n"; echo
