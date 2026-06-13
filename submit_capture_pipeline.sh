#!/bin/bash
#$ -cwd
#$ -o /u/scratch/a/aliceg/logs/misc/capture_pipeline.$JOB_ID.out
#$ -j y
# Capture foldseek's exact internal structure-search pipeline (the generated
# structuresearch.sh) so we can replicate prefilter+align+convert per shard.
#$ -l h_rt=0:30:00,h_data=16G,highp
#$ -pe shared 4

. /u/local/Modules/default/init/modules.sh
source ~/.bashrc
set -euo pipefail

FS="$HOME/prostT5-runner/foldseek/bin/foldseek"
SCR=/u/scratch/a/aliceg
WORK="$SCR/capture_pipeline"
rm -rf "$WORK"; mkdir -p "$WORK"

# Tiny query: first 5 entries of the existing chunk05 query DB.
awk 'NR<=5 {print $1}' "$SCR/rtx_vs_h200_test/rtx_chunk05/query_db/queryDB.index" > "$WORK/keys.txt"
"$FS" createsubdb "$WORK/keys.txt" "$SCR/rtx_vs_h200_test/rtx_chunk05/query_db/queryDB" "$WORK/tinyQ"

# CPU search vs one shard, keep tmp so we can read the generated pipeline script.
"$FS" search "$WORK/tinyQ" "$SCR/rtx_shards_pad/shard_00_pad" "$WORK/res" "$WORK/tmp" \
    --gpu 0 --threads 4 --remove-tmp-files 0 || true

echo "===STRUCTURESEARCH_SH==="
find "$WORK/tmp" -name 'structuresearch.sh' -exec cat {} \; 2>/dev/null
echo "===CMD_FILES==="
ls -R "$WORK/tmp" | head -60
echo "=== capture done ==="
