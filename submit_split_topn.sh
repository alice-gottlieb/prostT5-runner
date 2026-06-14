#!/bin/bash
# Launch the foldseek_topn_pfam arrays split across both GPU pools on DISJOINT
# chunk ranges: chunks 1..K on H200 (exact, serial -- 1 GPU at a time), chunks
# K+1..200 on RTX2080Ti (sharded ~98%, many in parallel). Both run concurrently
# and write to separate output bases. Run from ~/prostT5-runner on the Hoffman
# login node (this only issues qsubs; no compute here).
#
# Usage: ./submit_split_topn.sh [K] [RTX_TC]
#   K       last chunk on H200 (default 75); RTX gets K+1..200.
#   RTX_TC  cap on concurrent RTX array tasks via -tc (default 20).
set -euo pipefail

K="${1:-75}"
RTX_TC="${2:-20}"
TOTAL=200

echo "=== H200 (exact): chunks 1-$K ==="
qsub -t "1-$K" submit_topn_pfam_array_full.sh

echo "=== RTX (sharded ~98%): chunks $((K + 1))-$TOTAL  (-tc $RTX_TC) ==="
qsub -t "$((K + 1))-$TOTAL" -tc "$RTX_TC" submit_topn_pfam_array_rtx.sh

echo "Submitted. Poll with: qstat -u $USER"
