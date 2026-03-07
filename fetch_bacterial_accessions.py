#!/usr/bin/env python3
"""
Fetch bacterial reference genome accessions from NCBI Datasets API
and split them into chunks for parallel downloading via qsub.

Usage:
    python fetch_bacterial_accessions.py --num-genomes 100 --chunks 5
    python fetch_bacterial_accessions.py --num-genomes all --chunks 10
"""

import argparse
import json
import math
import os
import requests


def fetch_bacterial_reference_accessions(api_key, num_genomes):
    """Fetch bacterial reference genome accessions from NCBI Datasets API."""
    url = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/2/dataset_report"
    headers = {"api-key": api_key, "Accept": "application/json"}
    params = {
        "filters.reference_only": "true",
        "filters.assembly_source": "refseq",
        "page_size": 1000,
    }

    accessions = []
    page_token = None

    while True:
        if page_token:
            params["page_token"] = page_token

        resp = requests.get(url, headers=headers, params=params)
        resp.raise_for_status()
        data = resp.json()

        for report in data.get("reports", []):
            acc = report.get("accession")
            if acc:
                accessions.append(acc)
                if num_genomes != "all" and len(accessions) >= num_genomes:
                    return accessions

        page_token = data.get("next_page_token")
        if not page_token:
            break

        print(f"  Fetched {len(accessions)} accessions so far...")

    return accessions


def split_and_write(accessions, chunks, output_dir):
    """Split accessions into chunk files."""
    os.makedirs(output_dir, exist_ok=True)
    chunk_size = math.ceil(len(accessions) / chunks)

    for i in range(chunks):
        chunk = accessions[i * chunk_size : (i + 1) * chunk_size]
        if not chunk:
            break
        chunk_file = os.path.join(output_dir, f"chunk_{i + 1:02d}.txt")
        with open(chunk_file, "w") as f:
            f.write("\n".join(chunk) + "\n")
        print(f"  {chunk_file}: {len(chunk)} accessions")

    return min(chunks, math.ceil(len(accessions) / chunk_size))


def main():
    parser = argparse.ArgumentParser(description="Fetch bacterial reference genome accessions from NCBI")
    parser.add_argument("--num-genomes", default="100",
                        help="Number of genomes to fetch, or 'all' (default: 100)")
    parser.add_argument("--chunks", type=int, default=5,
                        help="Number of chunks to split into (default: 5)")
    parser.add_argument("--api-key-file", default="ncbi_key.txt",
                        help="Path to file containing NCBI API key")
    parser.add_argument("--output-dir", default="download_chunks",
                        help="Directory for chunk files (default: download_chunks)")
    args = parser.parse_args()

    with open(args.api_key_file) as f:
        api_key = f.read().strip()

    num = args.num_genomes if args.num_genomes == "all" else int(args.num_genomes)

    print(f"Fetching {'all' if num == 'all' else num} bacterial reference genome accessions...")
    accessions = fetch_bacterial_reference_accessions(api_key, num)
    print(f"Got {len(accessions)} accessions.")

    actual_chunks = split_and_write(accessions, args.chunks, args.output_dir)
    print(f"\nReady to submit: qsub -t 1-{actual_chunks} submit_download_genomes.sh")


if __name__ == "__main__":
    main()
