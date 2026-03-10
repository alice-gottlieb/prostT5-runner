#!/usr/bin/env python3
"""Collect all accessions from assemblies_downloaded in task_*/metadata.json files."""

import json
import os
import sys
from pathlib import Path

def main():
    scratch = os.environ.get("SCRATCH")
    if not scratch:
        sys.exit("Error: $SCRATCH environment variable is not set")

    base_dir = Path(scratch) / "all_3dis"
    if not base_dir.is_dir():
        sys.exit(f"Error: {base_dir} does not exist")

    from collections import Counter

    all_accessions = []

    # Search up to 3 levels deep for task_*/metadata.json
    for pattern in [
        "task_*/metadata.json",
        "*/task_*/metadata.json",
        "*/*/task_*/metadata.json",
    ]:
        for metadata_file in base_dir.glob(pattern):
            with open(metadata_file, "r") as f:
                data = json.load(f)
            all_accessions.extend(data.get("assemblies_downloaded", []))

    counts = Counter(all_accessions)
    unique = sorted(counts)
    duplicates = sorted(acc for acc, count in counts.items() if count > 1)

    with open("all_accessions.txt", "w") as f:
        for acc in unique:
            f.write(acc + "\n")

    with open("duplicate_accessions.txt", "w") as f:
        for acc in duplicates:
            f.write(f"{acc}\t{counts[acc]}\n")

    print(f"Collected {len(unique)} unique accessions into all_accessions.txt")
    print(f"Found {len(duplicates)} duplicate accessions in duplicate_accessions.txt")


if __name__ == "__main__":
    main()
