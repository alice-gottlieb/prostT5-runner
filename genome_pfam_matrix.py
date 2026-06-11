#!/usr/bin/env python3
"""Build a genome x pfam hit-count matrix from foldseek_topn_pfam.py's hits.tsv.

Each cell [genome, pfam] is the number of foldseek hits attributed to that
genome for that Pfam domain. The genome and Pfam both live in the `query`
column (col 1), which foldseek_topn_pfam.py builds as:

    {pfam}|{genome}|{target}|rank{rank}

so we split on "|" and take the Pfam (field 0) and genome (field 1).

Pipeline:
  1. Read hits.tsv and parse genome + pfam out of the query column.
  2. Intermediate long table: group by genome, then by pfam, and count the hits
     in each (genome, pfam) group.
  3. Pivot that long table into a wide matrix: rows = genome, cols = pfam,
     values = hit count (missing combinations filled with 0).

Usage:
  uv run --with polars python genome_pfam_matrix.py <out_dir>/hits.tsv
  uv run --with polars python genome_pfam_matrix.py hits.tsv -o genome_pfam.tsv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import polars as pl

QUERY_COL = "query"
GENOME_COL = "genome"
PFAM_COL = "pfam"


def load_hits(path: Path) -> pl.DataFrame:
    """Read hits.tsv and parse genome + pfam out of the query column.

    The query column looks like "{pfam}|{genome}|{target}|rank{rank}", so the
    Pfam is pipe-field 0 and the genome is pipe-field 1.
    """
    fields = pl.col(QUERY_COL).str.split("|")
    return pl.read_csv(path, separator="\t").select(
        fields.list.get(1).alias(GENOME_COL),
        fields.list.get(0).alias(PFAM_COL),
    )


def count_by_genome_and_pfam(hits: pl.DataFrame) -> pl.DataFrame:
    """Intermediate: long table grouped by genome, then by pfam, with hit counts.

    One row per (genome, pfam) combination that actually occurs.
    """
    return (
        hits.group_by([GENOME_COL, PFAM_COL])
        .agg(pl.len().alias("hits"))
        .sort([GENOME_COL, PFAM_COL])
    )


def pivot_to_matrix(long: pl.DataFrame) -> pl.DataFrame:
    """Wide matrix: rows = genome, cols = pfam, values = hit count (0-filled)."""
    matrix = long.pivot(values="hits", index=GENOME_COL, on=PFAM_COL)
    pfam_cols = sorted(c for c in matrix.columns if c != GENOME_COL)
    return (
        matrix.select([GENOME_COL, *pfam_cols])
        .with_columns(pl.col(pfam_cols).fill_null(0))
        .sort(GENOME_COL)
    )


def build_matrix(hits_tsv: Path) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Return (intermediate long table, final wide matrix) for a hits.tsv."""
    hits = load_hits(hits_tsv)
    long = count_by_genome_and_pfam(hits)
    matrix = pivot_to_matrix(long)
    return long, matrix


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("hits_tsv", type=Path,
                        help="hits.tsv produced by foldseek_topn_pfam.py")
    parser.add_argument("-o", "--output", type=Path, default=None,
                        help="Output TSV path (default: genome_pfam_counts.tsv "
                             "next to hits.tsv).")
    args = parser.parse_args()

    output = args.output or args.hits_tsv.parent / "genome_pfam_counts.tsv"

    long, matrix = build_matrix(args.hits_tsv)

    # Observe both tables as we go.
    with pl.Config(tbl_rows=20, tbl_cols=20):
        print("Intermediate: grouped by genome, then by pfam (long form):")
        print(long)
        print(f"\nFinal: genome x pfam hit-count matrix "
              f"({matrix.height} genomes x {len(matrix.columns) - 1} pfams):")
        print(matrix)

    matrix.write_csv(output, separator="\t")
    print(f"\nWrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
