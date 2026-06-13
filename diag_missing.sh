#!/bin/bash
#$ -cwd
#$ -o /u/scratch/a/aliceg/logs/misc/diag_missing.$JOB_ID.out
#$ -j y
#$ -l h_rt=0:30:00,h_data=16G,highp
#$ -pe shared 4

. /u/local/Modules/default/init/modules.sh
source ~/.bashrc
set -euo pipefail

FS="$HOME/prostT5-runner/foldseek/bin/foldseek"
SCR=/u/scratch/a/aliceg
W="$SCR/capture_pipeline/proto"
FULL="$SCR/all_3dis_fully_merged_2026-05-06/mergedDB"
Q="$SCR/capture_pipeline/tinyQ"
QID="PF12574.14|GCF_000008045.1|WP_011190939.1|rank2"

# Resolve the merged prefilter back to accessions so we can check membership.
"$FS" createtsv "$Q" "$FULL" "$W/pref_merged" "$W/pref_merged.tsv" --threads 4

awk -F'\t' -v q="$QID" '$1==q{print $2}' "$W/results2.tsv" | sort -u > "$W/new_targets.txt"
awk -F'\t' -v q="$QID" '$1==q{print $2}' "$SCR/rtx_vs_h200_test/h200_chunk05/results.full.tsv" | sort -u > "$W/ref_targets.txt"
awk -F'\t' -v q="$QID" '$1==q{print $2}' "$W/pref_merged.tsv" | sort -u > "$W/pref_targets.txt"

echo "new aligned targets:  $(wc -l < "$W/new_targets.txt")"
echo "ref aligned targets:  $(wc -l < "$W/ref_targets.txt")"
echo "pref_merged targets:  $(wc -l < "$W/pref_targets.txt")"

echo "=== targets in REF but not in NEW (the deficit) ==="
comm -23 "$W/ref_targets.txt" "$W/new_targets.txt" > "$W/missing.txt"
echo "count: $(wc -l < "$W/missing.txt")"
echo "--- are the missing targets in pref_merged? ---"
in_pref=$(comm -12 "$W/missing.txt" "$W/pref_targets.txt" | wc -l)
echo "missing targets present in pref_merged: $in_pref / $(wc -l < "$W/missing.txt")"
echo "--- missing targets with their H200 e-value/bits ---"
awk -F'\t' -v q="$QID" 'NR==FNR{m[$1];next} ($1==q && ($2 in m)){print $2"\t"$11"\t"$12}' \
    "$W/missing.txt" "$SCR/rtx_vs_h200_test/h200_chunk05/results.full.tsv" | sort -t$'\t' -k2 -g | head -20
echo "=== targets in NEW but not REF ==="
comm -13 "$W/ref_targets.txt" "$W/new_targets.txt" | head
echo "=== diag_missing done ==="
