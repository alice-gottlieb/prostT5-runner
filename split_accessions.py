"""
Split an accessions file into N roughly equal chunks for job array submission.

Usage:
    python split_accessions.py reference_genomes.txt --chunks 12
    python split_accessions.py reference_genomes.txt --chunks 12 --output-dir my_chunks
"""

import argparse
import os
import math


def main():
    parser = argparse.ArgumentParser(
        description="Split an accessions file into N chunks for parallel processing."
    )
    parser.add_argument("accession_file", help="Text file with one accession per line")
    parser.add_argument("--chunks", "-n", type=int, default=12, help="Number of chunks (default: 12)")
    parser.add_argument("--output-dir", "-o", default="chunks", help="Output directory (default: chunks)")
    args = parser.parse_args()

    # Read accessions (skip blanks and comments)
    with open(args.accession_file) as f:
        accessions = [line.strip() for line in f if line.strip() and not line.startswith("#")]

    n = len(accessions)
    k = min(args.chunks, n)
    chunk_size = math.ceil(n / k)

    os.makedirs(args.output_dir, exist_ok=True)

    for i in range(k):
        chunk = accessions[i * chunk_size : (i + 1) * chunk_size]
        chunk_file = os.path.join(args.output_dir, f"chunk_{i+1:02d}.txt")
        with open(chunk_file, "w") as f:
            f.write("\n".join(chunk) + "\n")
        print(f"  {chunk_file}: {len(chunk)} accessions")

    print(f"\nSplit {n} accessions into {k} chunks in {args.output_dir}/")


if __name__ == "__main__":
    main()
