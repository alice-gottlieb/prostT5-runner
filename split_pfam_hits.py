#!/usr/bin/env python3
"""Split a top-N Pfam-hits TSV into per-Pfam-subset chunk TSVs for array jobs.

foldseek_topn_pfam.py processes one input TSV at a time, building a query DB from
its 3Di codes. To fan that work out across an SGE GPU array job, we slice the input
by Pfam accession: each chunk is a smaller, self-contained input TSV (full header
plus every row whose Pfam falls in that chunk). Array task N then runs
foldseek_topn_pfam.py on pfam_chunk_NN.tsv.

Pfam accessions are read from query_pfam_accession (preferred) or domain_accession,
matching foldseek_topn_pfam.py. Distinct accessions are kept in first-seen order and
assigned to chunks round-robin so chunk sizes stay balanced.

Examples:
  # Test split: first 6 Pfams across 3 chunks (2 Pfams each).
  uv run python split_pfam_hits.py all_pfam_top50_hits.tsv \
      --num-chunks 3 --max-pfams 6 --out-dir pfam_chunks_test

  # Full split of every Pfam into 200 chunks.
  uv run python split_pfam_hits.py all_pfam_top50_hits.tsv \
      --num-chunks 200 --out-dir pfam_chunks

  # Explicit Pfam allow-list (overrides --max-pfams).
  uv run python split_pfam_hits.py all_pfam_top50_hits.tsv \
      --num-chunks 3 --pfams PF00001.21,PF00004.36,PF00069.30 --out-dir pfam_chunks_test
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path

PFAM_COLUMNS = ("query_pfam_accession", "domain_accession")


def log(message: str) -> None:
    print(message, flush=True)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("input_tsv", type=Path,
                   help="Top-N Pfam-hits TSV (e.g. all_pfam_top50_hits.tsv).")
    p.add_argument("--num-chunks", type=int, required=True,
                   help="Number of chunk TSVs to write (= number of array tasks).")
    p.add_argument("--out-dir", type=Path, required=True,
                   help="Directory to write pfam_chunk_NN.tsv files into.")
    p.add_argument("--max-pfams", type=int, default=None,
                   help="Keep only the first K distinct Pfam accessions (for small tests).")
    p.add_argument("--pfams", default=None,
                   help="Comma-separated Pfam allow-list. Overrides --max-pfams.")
    return p.parse_args()


def pick_pfam_column(fieldnames: list[str]) -> str:
    for column in PFAM_COLUMNS:
        if column in fieldnames:
            return column
    raise SystemExit(
        f"ERROR: input TSV has no Pfam column; expected one of {list(PFAM_COLUMNS)}."
    )


def choose_pfams(input_tsv: Path, pfam_column: str,
                 allow_list: set[str] | None, max_pfams: int | None) -> list[str]:
    """Distinct Pfam accessions in first-seen order, after applying any filters."""
    ordered: list[str] = []
    seen: set[str] = set()
    with input_tsv.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            pfam = row[pfam_column]
            if not pfam or pfam in seen:
                continue
            if allow_list is not None and pfam not in allow_list:
                continue
            seen.add(pfam)
            ordered.append(pfam)
            if allow_list is None and max_pfams is not None and len(ordered) >= max_pfams:
                break
    if not ordered:
        raise SystemExit("ERROR: no Pfam accessions matched the requested filters.")
    return ordered


def assign_chunks(pfams: list[str], num_chunks: int) -> dict[str, int]:
    """Round-robin each Pfam to a chunk index (1-based) for balanced chunk sizes."""
    return {pfam: (i % num_chunks) + 1 for i, pfam in enumerate(pfams)}


def write_chunks(input_tsv: Path, pfam_column: str, chunk_of_pfam: dict[str, int],
                 num_chunks: int, out_dir: Path) -> Counter:
    out_dir.mkdir(parents=True, exist_ok=True)
    chunk_paths = {n: out_dir / f"pfam_chunk_{n:02d}.tsv" for n in range(1, num_chunks + 1)}

    with input_tsv.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames or []

        files = {n: path.open("w", newline="") for n, path in chunk_paths.items()}
        writers = {n: csv.DictWriter(f, fieldnames=header, delimiter="\t",
                                     lineterminator="\n") for n, f in files.items()}
        for writer in writers.values():
            writer.writeheader()

        rows_per_chunk: Counter = Counter()
        try:
            for row in reader:
                chunk = chunk_of_pfam.get(row[pfam_column])
                if chunk is None:
                    continue
                writers[chunk].writerow(row)
                rows_per_chunk[chunk] += 1
        finally:
            for f in files.values():
                f.close()
    return rows_per_chunk


def main() -> int:
    args = parse_args()
    if args.num_chunks < 1:
        raise SystemExit("ERROR: --num-chunks must be at least 1.")

    allow_list = {p.strip() for p in args.pfams.split(",") if p.strip()} if args.pfams else None

    with args.input_tsv.open(newline="") as handle:
        pfam_column = pick_pfam_column(csv.DictReader(handle, delimiter="\t").fieldnames or [])
    log(f"Using Pfam column: {pfam_column}")

    pfams = choose_pfams(args.input_tsv, pfam_column, allow_list, args.max_pfams)
    log(f"Selected {len(pfams)} distinct Pfam accessions across {args.num_chunks} chunks.")

    chunk_of_pfam = assign_chunks(pfams, args.num_chunks)
    rows_per_chunk = write_chunks(args.input_tsv, pfam_column, chunk_of_pfam,
                                  args.num_chunks, args.out_dir)

    pfams_per_chunk: Counter = Counter(chunk_of_pfam.values())
    log(f"\nWrote chunks to {args.out_dir}:")
    for chunk in range(1, args.num_chunks + 1):
        members = [p for p in pfams if chunk_of_pfam[p] == chunk]
        log(f"  pfam_chunk_{chunk:02d}.tsv: "
            f"{pfams_per_chunk[chunk]} pfams, {rows_per_chunk[chunk]} rows "
            f"[{', '.join(members)}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
