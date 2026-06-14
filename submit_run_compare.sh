#!/bin/bash
#$ -cwd
#$ -o /u/scratch/a/aliceg/logs/misc/run_compare.$JOB_ID.out
#$ -j y
# Run the RTX-sharded vs full-DB equivalence comparison. The genome_counts
# matrices unpivot to tens of millions of cells, so this needs real host RAM --
# not a login-node job. Pass REF_DIR / TEST_DIR via `qsub -v`.
#$ -l h_rt=1:00:00,h_data=8G,highp
#$ -pe shared 8
#$ -M $USER@ucla.edu
#$ -m ea

. /u/local/Modules/default/init/modules.sh
source ~/.bashrc
set -euo pipefail

REF_DIR="${REF_DIR:?set REF_DIR}"
TEST_DIR="${TEST_DIR:?set TEST_DIR}"

uv run python -u ~/prostT5-runner/compare_rtx_h200.py "$REF_DIR" "$TEST_DIR"
echo "=== compare done ==="
