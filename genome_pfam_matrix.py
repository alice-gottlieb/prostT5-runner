#!/usr/bin/env python3
"""Build a pfam x genome hit-count matrix from foldseek_topn_pfam.py's hits.tsv.

Each cell [pfam, genome] is the number of foldseek hits for that Pfam domain
that landed in that *target* genome — i.e. how widely each Pfam's 3Di code is
found across the searched genomes. Genomes are the columns so each row is one
Pfam's profile across all genomes.

hits.tsv carries two ready-made columns we use directly:
  - `target_genome`:  the genome(s) the hit target belongs to. A multispecies
    target lives in several genomes, so foldseek_topn_pfam.py writes them
    ";"-joined; we split and count the hit toward each one.
  - `pfam_accession`: the Pfam the query came from.

(Earlier this script parsed the genome out of the `query` column, which is the
query's *source* genome, not the target genome it matched — collapsing every
hit onto the one genome the query came from.)

Pipeline:
  1. Read hits.tsv; take pfam_accession + target_genome, split/explode the
     ";"-joined target genomes, and drop unmapped ("NA") genomes.
  2. Intermediate long table: group by pfam, then by genome, and count the hits
     in each (pfam, genome) group.
  3. Pivot that long table into a wide matrix: rows = pfam, cols = genome,
     values = hit count (missing combinations filled with 0).

Usage:
  uv run --with polars python genome_pfam_matrix.py <out_dir>/hits.tsv
  uv run --with polars python genome_pfam_matrix.py hits.tsv -o genome_pfam.tsv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import polars as pl

TARGET_GENOME_COL = "target_genome"
PFAM_ACCESSION_COL = "pfam_accession"
GENOME_COL = "genome"
PFAM_COL = "pfam"


def load_hits(path: Path) -> pl.DataFrame:
    """Read hits.tsv and pull the Pfam + target genome(s) for each hit.

    `target_genome` may list several ";"-joined genomes for one multispecies
    target, so we split and explode into one row per (pfam, genome). Unmapped
    targets ("NA") are dropped.
    """
    return (
        pl.read_csv(path, separator="\t")
        .select(
            pl.col(TARGET_GENOME_COL).str.split(";").alias(GENOME_COL),
            pl.col(PFAM_ACCESSION_COL).alias(PFAM_COL),
        )
        .explode(GENOME_COL)
        .filter(pl.col(GENOME_COL) != "NA")
    )


def count_by_pfam_and_genome(hits: pl.DataFrame) -> pl.DataFrame:
    """Intermediate: long table grouped by pfam, then by genome, with hit counts.

    One row per (pfam, genome) combination that actually occurs.
    """
    return (
        hits.group_by([PFAM_COL, GENOME_COL])
        .agg(pl.len().alias("hits"))
        .sort([PFAM_COL, GENOME_COL])
    )


def pivot_to_matrix(long: pl.DataFrame) -> pl.DataFrame:
    """Wide matrix: rows = pfam, cols = genome, values = hit count (0-filled)."""
    matrix = long.pivot(values="hits", index=PFAM_COL, on=GENOME_COL)
    genome_cols = sorted(c for c in matrix.columns if c != PFAM_COL)
    return (
        matrix.select([PFAM_COL, *genome_cols])
        .with_columns(pl.col(genome_cols).fill_null(0))
        .sort(PFAM_COL)
    )


def build_matrix(hits_tsv: Path) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Return (intermediate long table, final wide matrix) for a hits.tsv."""
    hits = load_hits(hits_tsv)
    long = count_by_pfam_and_genome(hits)
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
        print("Intermediate: grouped by pfam, then by genome (long form):")
        print(long)
        print(f"\nFinal: pfam x genome hit-count matrix "
              f"({matrix.height} pfams x {len(matrix.columns) - 1} genomes):")
        print(matrix)

    matrix.write_csv(output, separator="\t")
    print(f"\nWrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
