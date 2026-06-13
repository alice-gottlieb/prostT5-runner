#!/bin/bash
#$ -cwd
#$ -o /u/scratch/a/aliceg/logs/misc/validate_remap.$JOB_ID.out
#$ -j y
# Validate the accession-merge + full-DB-keyspace rebuild on the 5 tiny queries:
# per-shard createtsv -> build_merged_prefilter (full-DB keys) -> structurealign
# vs FULL DB -> convertalis, then diff against the H200 reference rows.
#$ -l h_rt=1:00:00,h_data=24G,highp
#$ -pe shared 4

. /u/local/Modules/default/init/modules.sh
source ~/.bashrc
set -euo pipefail

FS="$HOME/prostT5-runner/foldseek/bin/foldseek"
SCR=/u/scratch/a/aliceg
SHARDS="$SCR/rtx_shards_pad"
FULL="$SCR/all_3dis_fully_merged_2026-05-06/mergedDB"
Q="$SCR/capture_pipeline/tinyQ"
W="$SCR/capture_pipeline/proto"
T=4
STRUALIGN_PAR="--tmscore-threshold 0 --tmscore-threshold-mode 0 --lddt-threshold 0 --sort-by-structure-bits 1 --alignment-type 2 --exact-tmscore 0 --sub-mat aa:3di.out,nucl:3di.out -a 0 --alignment-mode 3 --alignment-output-mode 0 --wrapped-scoring 0 -e 10 --min-seq-id 0 --min-aln-len 0 --seq-id-mode 0 --alt-ali 0 -c 0 --cov-mode 0 --max-seq-len 65535 --comp-bias-corr 1 --comp-bias-corr-scale 0.5 --max-rejected 2147483647 --max-accept 2147483647 --add-self-matches 0 --db-load-mode 0 --score-bias 0 --realign 0 --gap-open aa:10,nucl:10 --gap-extend aa:1,nucl:1 --zdrop 40 --threads $T -v 3"
COLS="query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"

echo "=== createtsv each shard prefilter -> accession TSVs ==="
TSVS=()
for k in 00 01 02 03 04 05; do
    "$FS" createtsv "$Q" "$SHARDS/shard_${k}_pad" "$W/pref_${k}" "$W/pref_${k}.tsv" --threads $T
    TSVS+=("$W/pref_${k}.tsv")
done

echo "=== build merged prefilter in full-DB keyspace ==="
uv run python -u ~/prostT5-runner/build_merged_prefilter.py \
    --query-db "$Q" --full-db "$FULL" --out "$W/pref_merged" --max-seqs 1000 "${TSVS[@]}"

echo "=== structurealign vs FULL DB ==="
# shellcheck disable=SC2086
"$FS" structurealign "$Q" "$FULL" "$W/pref_merged" "$W/aln2" $STRUALIGN_PAR
"$FS" convertalis "$Q" "$FULL" "$W/aln2" "$W/results2.tsv" --format-output "$COLS" --threads $T

echo "=== per-query hit counts: NEW sharded vs H200 ref ==="
cut -f1 "$W/results2.tsv" | sort -u > "$W/q2.txt"
for q in $(cat "$W/q2.txt"); do
    n_new=$(awk -F'\t' -v q="$q" '$1==q' "$W/results2.tsv" | wc -l)
    n_ref=$(awk -F'\t' -v q="$q" '$1==q' "$SCR/rtx_vs_h200_test/h200_chunk05/results.full.tsv" | wc -l)
    echo "  $q  new=$n_new  ref=$n_ref"
done

echo "=== exact (target,evalue,bits) agreement for one query ==="
Q1=$(head -1 "$W/q2.txt")
echo "query: $Q1"
echo "-- NEW (top 8 by bits) --"
awk -F'\t' -v q="$Q1" '$1==q{print $2"\t"$11"\t"$12}' "$W/results2.tsv" | sort -t$'\t' -k3 -nr | head -8
echo "-- H200 (top 8 by bits) --"
awk -F'\t' -v q="$Q1" '$1==q{print $2"\t"$11"\t"$12}' "$SCR/rtx_vs_h200_test/h200_chunk05/results.full.tsv" | sort -t$'\t' -k3 -nr | head -8
echo "=== validate done ==="
