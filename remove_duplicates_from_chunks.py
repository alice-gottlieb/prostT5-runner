#!/usr/bin/env python3
"""Remove completed accessions (from all_3di_chunks_400_split) from new_chunks_20."""

import sys
from pathlib import Path


def load_accessions_from_chunks(chunk_dir):
    """Load all accessions from chunk files in a directory."""
    accessions = set()
    for chunk_file in chunk_dir.iterdir():
        if chunk_file.is_file():
            with open(chunk_file) as f:
                for line in f:
                    line = line.strip()
                    if line:
                        accessions.add(line)
    return accessions


def main():
    completed_dir = Path("/u/home/a/aliceg/prostT5-runner/all_3di_chunks_400_split")
    if not completed_dir.is_dir():
        sys.exit(f"Error: {completed_dir} does not exist")

    target_dir = Path("/u/scratch/a/aliceg/new_chunks_20")
    if not target_dir.is_dir():
        sys.exit(f"Error: {target_dir} does not exist")

    # Load all accessions from the completed chunks directory
    completed = load_accessions_from_chunks(completed_dir)
    print(f"Loaded {len(completed)} accessions from {completed_dir}")

    total_removed = 0
    files_modified = 0

    for chunk_file in sorted(target_dir.iterdir()):
        if not chunk_file.is_file():
            continue

        with open(chunk_file) as f:
            lines = f.readlines()

        filtered = [line for line in lines if line.strip() not in completed]
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
