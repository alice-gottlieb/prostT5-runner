#!/usr/bin/env bash
# Run hmmsearch with Pfam-A.hmm against a single .faa file.
# One per-domain table (.domtblout) and one human-readable output file are written.
#
# Usage:
#   ./hmmsearch_single.sh <input_faa> <output_dir> <pfam_hmm> [cpus]
#
# Example:
#   ./hmmsearch_single.sh proteins/sample.faa pfam_hits/ pfam_data/Pfam-A.hmm 8

# set -euo pipefail

if [[ $# -lt 3 || $# -gt 4 ]]; then
    echo "Usage: $0 <input_faa> <output_dir> <pfam_hmm> [cpus]" >&2
    exit 1
fi

input_faa=$1
output_dir=$2
pfam_hmm=$3
cpus=${4:-4}

if ! [[ "$cpus" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: cpus must be a positive integer: $cpus" >&2
    exit 1
fi

if [[ ! -f "$input_faa" ]]; then
    echo "Error: input .faa file not found: $input_faa" >&2
    exit 1
fi

if [[ ! -f "$pfam_hmm" ]]; then
    echo "Error: HMM file not found: $pfam_hmm" >&2
    exit 1
fi

# hmmpress writes Pfam-A.hmm.h3{f,i,m,p} next to the .hmm file. Required once.
for ext in h3f h3i h3m h3p; do
    if [[ ! -f "$pfam_hmm.$ext" ]]; then
        echo "Missing $pfam_hmm.$ext - running hmmpress once..."
        hmmpress "$pfam_hmm"
        break
    fi
done

mkdir -p "$output_dir"

filename=$(basename "$input_faa")
base=${filename%.faa}
if [[ "$base" == "$filename" ]]; then
    base=${filename%.*}
fi

domtbl="$output_dir/$base.domtblout"
humanreadable="$output_dir/$base.humanreadable.out"

if [[ -s "$domtbl" && -s "$humanreadable" ]]; then
    echo "skip (exists): $base.domtblout and $base.humanreadable.out"
    echo "Done. Tables in $output_dir/"
    exit 0
fi

echo "hmmsearch: $filename -> $base.domtblout (cpus=$cpus)"
hmmsearch \
    --cut_ga \
    --cpu "$cpus" \
    --domtblout "$domtbl" \
    -o "$humanreadable" \
    "$pfam_hmm" \
    "$input_faa"

status=$?
if [[ "$status" -ne 0 ]]; then
    echo "Error: hmmsearch failed for $input_faa" >&2
    exit "$status"
fi

echo "Done. Tables in $output_dir/"
