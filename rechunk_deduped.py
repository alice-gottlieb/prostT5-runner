#!/usr/bin/env python3
"""Read all accessions from new_chunks_20 and rechunk them into final_chunk_20 with 20 per file."""

import math
import os
import sys
from pathlib import Path


def main():
    # scratch = os.environ.get("SCRATCH")
    # if not scratch:
    #     sys.exit("Error: $SCRATCH environment variable is not set")

    source_dir = Path("/u/scratch/a/aliceg/new_chunks_20")
    if not source_dir.is_dir():
        sys.exit(f"Error: {source_dir} does not exist")

    # Collect all accessions from all chunk files
    accessions = []
    seen = set()
    for chunk_file in sorted(source_dir.iterdir()):
        if not chunk_file.is_file():
            continue
        with open(chunk_file) as f:
            for line in f:
                acc = line.strip()
                if acc and acc not in seen:
                    accessions.append(acc)
                    seen.add(acc)

    print(f"Read {len(accessions)} unique accessions from {source_dir}")

    if not accessions:
        print("No accessions found. Nothing to do.")
        return

    output_dir = Path("/u/scratch/a/aliceg/final_chunk_20")
    output_dir.mkdir(parents=True, exist_ok=True)

    chunk_size = 20
    num_chunks = math.ceil(len(accessions) / chunk_size)

    for i in range(num_chunks):
        chunk = accessions[i * chunk_size : (i + 1) * chunk_size]
        chunk_file = output_dir / f"chunk_{i + 1:02d}.txt"
        with open(chunk_file, "w") as f:
            f.write("\n".join(chunk) + "\n")

    print(f"Wrote {num_chunks} chunk files to {output_dir}")
    print(f"Ready to submit: qsub -t 1-{num_chunks} submit_3di_array.sh")


if __name__ == "__main__":
    main()
