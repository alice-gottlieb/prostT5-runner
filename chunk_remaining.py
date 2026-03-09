#!/usr/bin/env python3
"""
Create chunk files for genomes that haven't been processed yet.

Reads all .faa files in the genomes directory, checks which accessions
have already been processed (by reading metadata.json from completed tasks),
and splits the remaining accessions into chunks.

Usage:
    python chunk_remaining.py --genomes-dir /path/to/ncbi_genomes \
                              --results-dir /path/to/all_3dis \
                              --chunk-size 50 \
                              --output-dir chunks
"""

import argparse
import glob
import json
import os
import math


def get_available_accessions(genomes_dir):
    """Get all accessions from .faa files in the genomes directory."""
    accessions = set()
    for faa in glob.glob(os.path.join(genomes_dir, "*.faa")):
        acc = os.path.splitext(os.path.basename(faa))[0]
        accessions.add(acc)
    return accessions


def get_completed_accessions(results_dir):
    """Get accessions already processed from metadata.json in task_* dirs."""
    completed = set()
    for meta_path in glob.glob(os.path.join(results_dir, "task_*/metadata.json")):
        try:
            with open(meta_path) as f:
                meta = json.load(f)
            for acc in meta.get("assemblies_downloaded", []):
                completed.add(acc)
        except (json.JSONDecodeError, IOError) as e:
            print(f"  [WARN] Skipping {meta_path}: {e}")
    return completed


def main():
    parser = argparse.ArgumentParser(description="Chunk remaining unprocessed genomes")
    parser.add_argument("--genomes-dir", required=True,
                        help="Directory containing .faa files")
    parser.add_argument("--results-dir", required=True,
                        help="Directory containing task_*/metadata.json results")
    parser.add_argument("--chunk-size", type=int, default=50,
                        help="Number of accessions per chunk (default: 50)")
    parser.add_argument("--output-dir", default="chunks",
                        help="Directory for chunk files (default: chunks)")
    args = parser.parse_args()

    available = get_available_accessions(args.genomes_dir)
    print(f"Available genomes: {len(available)}")

    completed = get_completed_accessions(args.results_dir)
    print(f"Already processed: {len(completed)}")

    remaining = sorted(available - completed)
    print(f"Remaining: {len(remaining)}")

    if not remaining:
        print("Nothing left to process!")
        return

    os.makedirs(args.output_dir, exist_ok=True)
    num_chunks = math.ceil(len(remaining) / args.chunk_size)

    for i in range(num_chunks):
        chunk = remaining[i * args.chunk_size : (i + 1) * args.chunk_size]
        chunk_file = os.path.join(args.output_dir, f"chunk_{i + 1:02d}.txt")
        with open(chunk_file, "w") as f:
            f.write("\n".join(chunk) + "\n")
        print(f"  {chunk_file}: {len(chunk)} accessions")

    print(f"\nReady to submit: qsub -t 1-{num_chunks} submit_3di_array.sh")


if __name__ == "__main__":
    main()
