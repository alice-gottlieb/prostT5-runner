#!/bin/bash
#$ -cwd
#$ -o logs/makepadded.$JOB_ID.out
#$ -j y
#$ -l h_data=16G,h_rt=4:00:00
#$ -pe shared 8
#
# One-time prep: build a padded GPU target DB from the merged 3Di DB.
#
# foldseek `search --gpu 1` rejects an ordinary target DB ("not a valid GPU
# database") and asks for a makepaddedseqdb copy. This job builds that copy once;
# the array job (submit_topn_pfam_array.sh) then searches against PADDED_DB.

# --- Environment (source before `set -u`; /etc/bashrc uses an unbound PS1) ---
. /u/local/Modules/default/init/modules.sh
source ~/.bashrc
ulimit -c 0

set -euo pipefail

SOURCE_DB="$SCRATCH/all_3dis_fully_merged_2026-05-06/mergedDB"
PADDED_DB="$SCRATCH/all_3dis_fully_merged_2026-05-06_gpu_pad/mergedDB_pad"
FOLDSEEK_BIN="$HOME/prostT5-runner/foldseek/bin/foldseek"

mkdir -p "$(dirname "$PADDED_DB")" logs

echo "=== makepaddedseqdb: $SOURCE_DB -> $PADDED_DB ==="
"$FOLDSEEK_BIN" makepaddedseqdb "$SOURCE_DB" "$PADDED_DB" --threads 8

echo "=== padded DB files ==="
ls -la "$(dirname "$PADDED_DB")"
echo "=== done ==="
