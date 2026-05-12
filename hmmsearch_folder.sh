#!/usr/bin/env bash
# Run hmmsearch with Pfam-A.hmm against every .faa file in an input folder.
# One per-domain table (.domtblout) is written per input file.
#
# Usage:
#   ./hmmsearch_folder.sh <input_dir> <output_dir> <pfam_hmm> [cpus_per_job] [max_jobs]
#
# Example:
#   ./hmmsearch_folder.sh proteins/ pfam_hits/ pfam_data/Pfam-A.hmm 8 2

# set -euo pipefail

if [[ $# -lt 3 || $# -gt 5 ]]; then
    echo "Usage: $0 <input_dir> <output_dir> <pfam_hmm> [cpus_per_job] [max_jobs]" >&2
    exit 1
fi

input_dir=$1
output_dir=$2
pfam_hmm=$3
cpus_per_job=${4:-4}
max_jobs=${5:-1}

if ! [[ "$cpus_per_job" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: cpus_per_job must be a positive integer: $cpus_per_job" >&2
    exit 1
fi

if ! [[ "$max_jobs" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: max_jobs must be a positive integer: $max_jobs" >&2
    exit 1
fi

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

echo "Found ${#faa_files[@]} .faa files. Output: $output_dir (cpus_per_job=$cpus_per_job, max_jobs=$max_jobs)"

pids=()
for faa in "${faa_files[@]}"; do
    base=$(basename "$faa" .faa)
    domtbl="$output_dir/$base.domtblout"
    humanreadable="$output_dir/$base.humanreadable.out"

    if [[ -s "$domtbl" && -s "$humanreadable" ]]; then
        echo "  skip (exists): $base.domtblout and $base.humanreadable.out"
        continue
    fi

    echo "  hmmsearch: $base.faa -> $base.domtblout"
    hmmsearch \
        --cut_ga \
        --cpu "$cpus_per_job" \
        --domtblout "$domtbl" \
        -o "$humanreadable" \
        "$pfam_hmm" \
        "$faa" &
    pids+=("$!")

    while [[ $(jobs -rp | wc -l) -ge "$max_jobs" ]]; do
        sleep 1
    done
done

status=0
for pid in "${pids[@]}"; do
    if ! wait "$pid"; then
        status=1
    fi
done

if [[ "$status" -ne 0 ]]; then
    echo "Error: one or more hmmsearch jobs failed" >&2
    exit "$status"
fi

echo "Done. Tables in $output_dir/"
