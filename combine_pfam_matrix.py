#!/usr/bin/env python3
"""Stitch per-chunk Pfam x genome matrices from one or more output bases into one
full matrix, recording which base (GPU pool) produced each chunk.

The chunks partition the Pfams disjointly, so each chunk's genome_pfam_matrix.tsv
(rows = pfam, cols = genome) holds a distinct set of Pfam rows. The full matrix is
the row-wise stack of every chunk, aligning genome columns across chunks (union;
genomes absent from a chunk are 0). This lets an H200 run (exact) and an RTX run
(sharded, ~98%) cover different chunks and still combine into one matrix.

Usage:
  uv run python combine_pfam_matrix.py \
      --base H200=$SCRATCH/foldseek_topn_pfam_array_full \
      --base RTX=$SCRATCH/foldseek_topn_pfam_array_rtx \
      --expect 200 -o $SCRATCH/full_matrix.tsv --source-out $SCRATCH/chunk_source.tsv
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import polars as pl

PFAM_COL = "pfam"
CHUNK_RE = re.compile(r"chunk_(\d+)")
MATRIX_NAME = "genome_pfam_matrix.tsv"


def discover_chunks(bases: list[tuple[str, Path]]) -> dict[str, tuple[str, Path]]:
    """chunk_id -> (label, matrix_path), erroring if a chunk shows up twice."""
    found: dict[str, tuple[str, Path]] = {}
    for label, base in bases:
        for matrix in sorted(base.glob(f"chunk_*/{MATRIX_NAME}")):
            m = CHUNK_RE.search(matrix.parent.name)
            if not m:
                continue
            chunk = m.group(1)
            if chunk in found:
                sys.exit(f"ERROR: chunk {chunk} found in both "
                         f"{found[chunk][0]} and {label} ({matrix}). Ranges must "
                         f"be disjoint.")
            found[chunk] = (label, matrix)
    return found


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", action="append", required=True, metavar="LABEL=DIR",
                    help="Output base and its label, e.g. H200=/path/array_full. "
                         "Repeat for each GPU pool.")
    ap.add_argument("--expect", type=int, default=200,
                    help="Expected number of chunks (1..N). Default 200.")
    ap.add_argument("-o", "--output", type=Path, required=True,
                    help="Full Pfam x genome matrix TSV.")
    ap.add_argument("--source-out", type=Path, default=None,
                    help="Provenance TSV: chunk <tab> source label.")
    ap.add_argument("--allow-incomplete", action="store_true",
                    help="Warn instead of error on missing chunks (partial runs).")
    args = ap.parse_args()

    bases = []
    for spec in args.base:
        label, _, path = spec.partition("=")
        if not path:
            sys.exit(f"ERROR: --base must be LABEL=DIR, got {spec!r}")
        bases.append((label, Path(path)))

    found = discover_chunks(bases)
    print(f"Discovered {len(found)} chunk matrices across {len(bases)} bases.")

    expected = {f"{i:02d}" for i in range(1, args.expect + 1)}
    missing = sorted(expected - set(found))
    if missing:
        msg = f"missing {len(missing)} chunks: {', '.join(missing)}"
        if args.allow_incomplete:
            print(f"WARNING: {msg}")
        else:
            sys.exit(f"ERROR: {msg} (use --allow-incomplete to combine anyway)")

    # Stack every chunk's matrix, aligning genome columns by name (union, 0-fill).
    frames = []
    for chunk in sorted(found):
        label, matrix = found[chunk]
        frames.append(pl.read_csv(matrix, separator="\t"))
    print(f"Concatenating {len(frames)} chunk matrices (diagonal join)...")
    full = pl.concat(frames, how="diagonal")

    genome_cols = sorted(c for c in full.columns if c != PFAM_COL)
    full = (full.select([PFAM_COL, *genome_cols])
            .with_columns(pl.col(genome_cols).fill_null(0))
            .sort(PFAM_COL))
    full.write_csv(args.output, separator="\t")
    print(f"Wrote {args.output} "
          f"({full.height:,} pfams x {len(genome_cols):,} genomes).")

    if args.source_out:
        src = pl.DataFrame({
            "chunk": [f"chunk_{c}" for c in sorted(found)],
            "source": [found[c][0] for c in sorted(found)],
        })
        src.write_csv(args.source_out, separator="\t")
        counts = dict(src["source"].value_counts().iter_rows())
        print(f"Wrote {args.source_out}: "
              + ", ".join(f"{k}={v}" for k, v in counts.items()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
