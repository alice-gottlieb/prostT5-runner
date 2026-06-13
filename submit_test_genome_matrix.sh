#!/bin/bash
# Test the genome_pfam_matrix.py step that the array scripts now run per chunk,
# using an already-produced hits.tsv (H200 pilot chunk_01) so we validate the
# exact command + cluster env without re-running the GPU foldseek search.
#$ -cwd
#$ -o logs/test_genome_matrix.$JOB_ID.out
#$ -j y
#$ -l h_rt=1:00:00,h_data=48G
#$ -M $USER@ucla.edu
#$ -m ea

. /u/local/Modules/default/init/modules.sh
source ~/.bashrc
set -euo pipefail

HITS="$SCRATCH/foldseek_topn_pfam_array_full/chunk_01/hits.tsv"
OUT="$SCRATCH/foldseek_topn_pfam_array_full/chunk_01/genome_pfam_matrix.tsv"

echo "=== Building Pfam x genome matrix from $HITS ==="
uv run python -u ~/prostT5-runner/genome_pfam_matrix.py "$HITS" -o "$OUT"

echo "=== Result ==="
echo "rows (pfams + header): $(wc -l < "$OUT")"
echo "cols (genomes + 1):    $(head -1 "$OUT" | tr '\t' '\n' | wc -l)"
echo "--- top-left preview ---"
cut -f1-4 "$OUT" | head -5
