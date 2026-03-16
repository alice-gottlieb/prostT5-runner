#!/usr/bin/env python3
"""Remove completed accessions (found in metadata.json files) from new_chunks_20 chunk files."""

import argparse
import json
import os
import sys
from pathlib import Path


def load_completed_accessions(base_dir):
    """Load all accessions from assemblies_downloaded in task_*/metadata.json files."""
    completed = set()
    for pattern in [
        "task_*/metadata.json",
        "*/task_*/metadata.json",
        "*/*/task_*/metadata.json",
    ]:
        for metadata_file in base_dir.glob(pattern):
            with open(metadata_file) as f:
                data = json.load(f)
            completed.update(data.get("assemblies_downloaded", []))
    return completed


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true", help="Print what would be removed without modifying files")
    args = parser.parse_args()

    scratch = os.environ.get("SCRATCH")
    if not scratch:
        sys.exit("Error: $SCRATCH environment variable is not set")

    base_dir = Path(scratch) / "all_3dis"
    if not base_dir.is_dir():
        sys.exit(f"Error: {base_dir} does not exist")

    completed = load_completed_accessions(base_dir)
    print(f"Loaded {len(completed)} completed accessions from metadata.json files")

    target_dir = Path("/u/scratch/a/aliceg/new_chunks_20")
    if not target_dir.is_dir():
        sys.exit(f"Error: {target_dir} does not exist")

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
            if not args.dry_run:
                if filtered:
                    with open(chunk_file, "w") as f:
                        f.writelines(filtered)
                else:
                    chunk_file.unlink()
            total_removed += removed
            files_modified += 1
            action = "would remove" if args.dry_run else ("deleted" if not filtered else "removed")
            print(f"  {chunk_file.name}: {action} {removed} accessions{' (file deleted)' if not args.dry_run and not filtered else ''}")

    if args.dry_run:
        print(f"\nDry run. Would remove {total_removed} lines from {files_modified} chunk files.")
    else:
        print(f"\nDone. Removed {total_removed} lines from {files_modified} chunk files.")


if __name__ == "__main__":
    main()
