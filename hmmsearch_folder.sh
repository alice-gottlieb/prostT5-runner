#!/usr/bin/env bash
# Run hmmsearch with Pfam-A.hmm against every .faa file in an input folder.
# One per-domain table (.domtblout) is written per input file.
#
# Usage:
#   ./run_hmmsearch_folder.sh <input_dir> <output_dir> <pfam_hmm> [cpus]
#
# Example:
#   ./run_hmmsearch_folder.sh proteins/ pfam_hits/ pfam_data/Pfam-A.hmm 8

set -euo pipefail

if [[ $# -lt 3 || $# -gt 4 ]]; then
    echo "Usage: $0 <input_dir> <output_dir> <pfam_hmm> [cpus]" >&2
    exit 1
fi

input_dir=$1
output_dir=$2
pfam_hmm=$3
cpus=${4:-4}

if [[ ! -d "$input_dir" ]]; then
    echo "Error: input dir not found: $input_dir" >&2
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

shopt -s nullglob
faa_files=("$input_dir"/*.faa)
if [[ ${#faa_files[@]} -eq 0 ]]; then
    echo "No .faa files found in $input_dir" >&2
    exit 1
fi

echo "Found ${#faa_files[@]} .faa files. Output: $output_dir (cpus=$cpus)"

for faa in "${faa_files[@]}"; do
    base=$(basename "$faa" .faa)
    domtbl="$output_dir/$base.domtblout"

    if [[ -s "$domtbl" ]]; then
        echo "  skip (exists): $base.domtblout"
        continue
    fi

    echo "  hmmsearch: $base.faa -> $base.domtblout"
    hmmsearch \
        --cut_ga \
        --cpu "$cpus" \
        --domtblout "$domtbl" \
        -o "$output_dir/$base.humanreadable.out" \
        "$pfam_hmm" \
        "$faa"
done

echo "Done. Tables in $output_dir/"
