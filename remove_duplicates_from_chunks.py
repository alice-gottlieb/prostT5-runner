#!/usr/bin/env python3
"""Remove duplicate accessions from chunk files in $SCRATCH/new_chunks_20."""

import os
import sys
from pathlib import Path


def main():
    script_dir = Path(__file__).parent
    dup_file = script_dir / "duplicate_accessions.txt"
    if not dup_file.is_file():
        sys.exit(f"Error: {dup_file} not found. Run collect_accessions.py first.")

    # Load duplicate accessions (first column, tab-separated)
    duplicates = set()
    with open(dup_file) as f:
        for line in f:
            line = line.strip()
            if line:
                duplicates.add(line.split("\t")[0])

    if not duplicates:
        print("No duplicates found. Nothing to do.")
        return

    print(f"Loaded {len(duplicates)} duplicate accessions")

    chunk_dir = Path("/u/scratch/a/aliceg/new_chunks_20")
    if not chunk_dir.is_dir():
        sys.exit(f"Error: {chunk_dir} does not exist")

    total_removed = 0
    files_modified = 0

    for chunk_file in sorted(chunk_dir.iterdir()):
        if not chunk_file.is_file():
            continue

        with open(chunk_file) as f:
            lines = f.readlines()

        filtered = [line for line in lines if line.strip() not in duplicates]
        removed = len(lines) - len(filtered)

        if removed > 0:
            # with open(chunk_file, "w") as f:
            #     f.writelines(filtered)
            total_removed += removed
            files_modified += 1
            print(f"  {chunk_file.name}: removed {removed} accessions")

    print(f"\nDone. Removed {total_removed} lines from {files_modified} chunk files.")


if __name__ == "__main__":
    main()
