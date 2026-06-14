#!/bin/bash
#$ -cwd
#$ -o /u/scratch/a/aliceg/logs/misc/combine_pfam_matrix.$JOB_ID.out
#$ -j y
# Combine per-chunk Pfam x genome matrices from the H200 (exact) and RTX (~98%)
# array output bases into one full matrix + a provenance manifest. The union
# matrix is large (~17k Pfams x tens of thousands of genomes), so ask for plenty
# of RAM. CPU-only -> highp is fine (never combine highp with a GPU request).
#$ -l h_rt=2:00:00,h_data=16G,highp
#$ -pe shared 8
#$ -M $USER@ucla.edu
#$ -m ea

. /u/local/Modules/default/init/modules.sh
source ~/.bashrc
set -euo pipefail

H200_BASE="$SCRATCH/foldseek_topn_pfam_array_full"
RTX_BASE="$SCRATCH/foldseek_topn_pfam_array_rtx"
OUT="$SCRATCH/foldseek_topn_pfam_full_matrix.tsv"
SRC="$SCRATCH/foldseek_topn_pfam_chunk_source.tsv"
EXPECT="${EXPECT:-200}"

uv run python -u ~/prostT5-runner/combine_pfam_matrix.py \
    --base "H200=$H200_BASE" \
    --base "RTX=$RTX_BASE" \
    --expect "$EXPECT" \
    -o "$OUT" \
    --source-out "$SRC"
echo "=== combine done ==="
