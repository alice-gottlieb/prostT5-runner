#!/bin/bash
#$ -cwd
#$ -o /u/scratch/a/aliceg/logs/misc/proto_pfmerge.$JOB_ID.out
#$ -j y
# Prototype the foldseek-native sharded pipeline on 5 tiny queries and check it
# reproduces the H200 full-DB e-values/bits. Per shard: GPU ungappedprefilter ->
# mergedbs (union) -> filterdb (global top-1000 by ungapped score) -> CPU
# structurealign vs the FULL DB (true full-DB e-values) -> convertalis.
#$ -l gpu,RTX2080Ti,cuda=1,h_rt=1:00:00,h_data=24G
#$ -M $USER@ucla.edu
#$ -m bea

. /u/local/Modules/default/init/modules.sh
source ~/.bashrc
ulimit -c 0
set -euo pipefail

FS="$HOME/prostT5-runner/foldseek/bin/foldseek"
SCR=/u/scratch/a/aliceg
SHARDS="$SCR/rtx_shards_pad"
FULL="$SCR/all_3dis_fully_merged_2026-05-06/mergedDB"
Q="$SCR/capture_pipeline/tinyQ"
W="$SCR/capture_pipeline/proto"; rm -rf "$W"; mkdir -p "$W"
T=4

# Params copied verbatim from foldseek search's own GPU pipeline (captured from
# the ungappedprefilter call and structuresearch.sh's structurealign call).
UNGAPPED_PAR="--sub-mat aa:3di.out,nucl:3di.out -c 0 -e 1.79769e+308 --cov-mode 0 --comp-bias-corr 1 --comp-bias-corr-scale 0.15 --min-ungapped-score 30 --max-seqs 1000 --db-load-mode 0 --gpu 1 --threads $T -v 3"
STRUALIGN_PAR="--tmscore-threshold 0 --tmscore-threshold-mode 0 --lddt-threshold 0 --sort-by-structure-bits 1 --alignment-type 2 --exact-tmscore 0 --sub-mat aa:3di.out,nucl:3di.out -a 0 --alignment-mode 3 --alignment-output-mode 0 --wrapped-scoring 0 -e 10 --min-seq-id 0 --min-aln-len 0 --seq-id-mode 0 --alt-ali 0 -c 0 --cov-mode 0 --max-seq-len 65535 --comp-bias-corr 1 --comp-bias-corr-scale 0.5 --max-rejected 2147483647 --max-accept 2147483647 --add-self-matches 0 --db-load-mode 0 --score-bias 0 --realign 0 --gap-open aa:10,nucl:10 --gap-extend aa:1,nucl:1 --zdrop 40 --threads $T -v 3"
COLS="query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"

echo "=== per-shard GPU ungappedprefilter ==="
PREFS=()
for k in 00 01 02 03 04 05; do
    # shellcheck disable=SC2086
    "$FS" ungappedprefilter "${Q}_ss" "${SHARDS}/shard_${k}_pad_ss" "$W/pref_${k}" $UNGAPPED_PAR
    PREFS+=("$W/pref_${k}")
done

echo "=== mergedbs (union per query) ==="
"$FS" mergedbs "$Q" "$W/pref_all" "${PREFS[@]}"

echo "=== filterdb: global top-1000 by ungapped score (col 2) ==="
"$FS" filterdb "$W/pref_all" "$W/pref_top" --filter-column 2 --sort-entries 2 --extract-lines 1000 --threads $T

echo "=== structurealign vs FULL DB (true full-DB e-values) ==="
# shellcheck disable=SC2086
"$FS" structurealign "$Q" "$FULL" "$W/pref_top" "$W/aln" $STRUALIGN_PAR

echo "=== convertalis ==="
"$FS" convertalis "$Q" "$FULL" "$W/aln" "$W/proto_results.tsv" --format-output "$COLS" --threads $T

echo "=== PROTO results (query,target,evalue,bits) ==="
sort "$W/proto_results.tsv" | awk -F'\t' '{print $1"\t"$2"\t"$11"\t"$12}' | head -20
echo "=== H200 reference rows for the SAME queries ==="
cut -f1 "$W/proto_results.tsv" | sort -u > "$W/proto_qids.txt"
awk -F'\t' 'NR==FNR{q[$1];next} ($1 in q){print $1"\t"$2"\t"$11"\t"$12}' \
    "$W/proto_qids.txt" "$SCR/rtx_vs_h200_test/h200_chunk05/results.full.tsv" | sort | head -20
echo "=== proto done ==="
