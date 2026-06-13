#!/bin/bash
#$ -cwd
#$ -o /u/scratch/a/aliceg/logs/misc/build_rtx_shards.$JOB_ID.out
#$ -j y
# Build NUM_SHARDS disjoint padded target DBs that together partition the full
# 90.1M-sequence merged 3Di DB. Each shard's padded _ss must fit an RTX2080Ti's
# ~10.5G VRAM: the full _ss is 29.4G and OOMs the GPU prefilter (--split does NOT
# help -- the padded target is loaded whole), but the proven 1/6 subset (~4.9G
# _ss) fits. 6 disjoint shards reproduce that subset size while covering ALL
# targets, so a search over all shards (merged) covers the whole DB.
# CPU-only build (createsubdb + makepaddedseqdb); never combine highp with a GPU.
#$ -l h_rt=8:00:00,h_data=6G,highp
#$ -pe shared 8
#$ -M $USER@ucla.edu
#$ -m ea

. /u/local/Modules/default/init/modules.sh
source ~/.bashrc
set -euo pipefail

FS="$HOME/prostT5-runner/foldseek/bin/foldseek"
SRC="$SCRATCH/all_3dis_fully_merged_2026-05-06/mergedDB"
OUTDIR="$SCRATCH/rtx_shards_pad"
NUM_SHARDS=6
mkdir -p "$OUTDIR"
cp "$(realpath "$0")" "$OUTDIR/$(basename "$0")"

for k in $(seq 0 $((NUM_SHARDS - 1))); do
    kk=$(printf "%02d" "$k")
    echo "=== Building shard $kk / $((NUM_SHARDS - 1)) ==="

    # Disjoint partition: row (NR-1) goes to shard (NR-1) % NUM_SHARDS.
    KEYS="$OUTDIR/shard_${kk}_keys.txt"
    awk -v k="$k" -v n="$NUM_SHARDS" '((NR-1) % n) == k {print $1}' "${SRC}.index" > "$KEYS"
    echo "shard $kk keys: $(wc -l < "$KEYS")"

    # Subset the whole structure DB (sequence + 3di + headers), then pad for GPU.
    "$FS" createsubdb "$KEYS" "$SRC" "$OUTDIR/shard_${kk}" --id-mode 0
    "$FS" makepaddedseqdb "$OUTDIR/shard_${kk}" "$OUTDIR/shard_${kk}_pad" --threads 8
done

echo "=== padded shard _ss sizes (each must fit RTX VRAM ~10.5G) ==="
ls -lh "$OUTDIR"/shard_*_pad_ss
echo "=== done ==="
