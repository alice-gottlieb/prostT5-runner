#!/bin/bash
#$ -cwd
#$ -o /u/scratch/a/aliceg/logs/misc/debug_keys.$JOB_ID.out
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

echo "=== pref_top resolved via FULL DB headers (query, target, score) ==="
"$FS" createtsv "$Q" "$FULL" "$W/pref_top" "$W/pref_top_full.tsv" --threads 4
echo "rows: $(wc -l < "$W/pref_top_full.tsv")"
echo "--- top entries for first query ---"
head -15 "$W/pref_top_full.tsv"

echo ""
echo "=== pref_top resolved via SHARD_00 headers (same key -> which accession?) ==="
"$FS" createtsv "$Q" "$SCR/rtx_shards_pad/shard_00_pad" "$W/pref_top_shard.tsv" --threads 4 2>/dev/null || echo "(shard createtsv failed)"
head -5 "$W/pref_top_shard.tsv" 2>/dev/null

echo ""
echo "=== Does pref_top (via full DB) contain the known strong self-hit WP_011190939.1? ==="
grep -c "WP_011190939.1" "$W/pref_top_full.tsv" || echo "NOT FOUND in pref_top"
echo "=== Raw pref_top numeric keys for query 0 (first entry block) ==="
head -5 "$W/pref_top.0" 2>/dev/null | cat -A
echo "=== mergedDB.lookup: what accession is key 14772948? (sample from capture) ==="
awk -F'\t' '$1==14772948 {print; exit}' "$FULL.lookup" 2>/dev/null || echo "key not in lookup"
echo "=== shard_00_pad.lookup: what accession is key 14772948 there? ==="
awk -F'\t' '$1==14772948 {print; exit}' "$SCR/rtx_shards_pad/shard_00_pad.lookup" 2>/dev/null || echo "key not in shard lookup"
echo "=== debug done ==="
