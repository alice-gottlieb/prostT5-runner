#!/usr/bin/env python3
"""Compare a sharded-RTX foldseek_topn_pfam run against an H200 reference run.

Checks that the two runs produce statistically equivalent output:
  * hits.tsv          -> overlap of (query, target) hit pairs (Jaccard).
  * genome_counts.tsv -> cell-by-cell agreement of the genome x query matrix.

Usage:
  uv run python compare_rtx_h200.py REF_DIR TEST_DIR

REF_DIR and TEST_DIR each contain hits.tsv and genome_counts.tsv.
"""

from __future__ import annotations

import sys
from pathlib import Path

import polars as pl


def hit_pairs(hits_tsv: Path) -> pl.DataFrame:
    """Unique (query, target) pairs from a hits.tsv."""
    return (
        pl.scan_csv(hits_tsv, separator="\t")
        .select(["query", "target"])
        .unique()
        .collect()
    )


def compare_hits(ref: Path, test: Path) -> None:
    print("\n=== hits.tsv : (query, target) pair overlap ===")
    r = hit_pairs(ref / "hits.tsv")
    t = hit_pairs(test / "hits.tsv")
    rn, tn = r.height, t.height
    inter = r.join(t, on=["query", "target"], how="inner").height
    union = rn + tn - inter
    jac = inter / union if union else 1.0
    print(f"  ref pairs:        {rn:,}")
    print(f"  test pairs:       {tn:,}")
    print(f"  shared pairs:     {inter:,}")
    print(f"  ref-only pairs:   {rn - inter:,}")
    print(f"  test-only pairs:  {tn - inter:,}")
    print(f"  Jaccard overlap:  {jac:.6f}")


def load_counts(counts_tsv: Path) -> pl.DataFrame:
    """genome_counts.tsv as a long (genome, query, count) frame."""
    wide = pl.read_csv(counts_tsv, separator="\t")
    id_col = wide.columns[0]
    return wide.unpivot(index=id_col, variable_name="query", value_name="count") \
               .rename({id_col: "genome"})


def compare_counts(ref: Path, test: Path) -> None:
    print("\n=== genome_counts.tsv : cell-by-cell agreement ===")
    r = load_counts(ref / "genome_counts.tsv")
    t = load_counts(test / "genome_counts.tsv")

    # Full outer join on (genome, query); missing cells are 0.
    joined = (
        r.join(t, on=["genome", "query"], how="full", suffix="_test")
        .with_columns(
            pl.col("count").fill_null(0),
            pl.col("count_test").fill_null(0),
        )
    )
    # Drop cells that are 0 in both runs (the matrix is sparse; only compare
    # cells that at least one run populated).
    nz = joined.filter((pl.col("count") != 0) | (pl.col("count_test") != 0))
    total = nz.height
    diff = nz.with_columns((pl.col("count") - pl.col("count_test")).abs().alias("absdiff"))
    identical = diff.filter(pl.col("absdiff") == 0).height
    nonzero_diff = total - identical
    max_abs = diff.select(pl.col("absdiff").max()).item() if total else 0
    sum_abs = diff.select(pl.col("absdiff").sum()).item() if total else 0
    frac = identical / total if total else 1.0

    print(f"  non-empty cells (either run): {total:,}")
    print(f"  identical cells:              {identical:,}")
    print(f"  differing cells:              {nonzero_diff:,}")
    print(f"  fraction identical:           {frac:.6f}")
    print(f"  max abs difference:           {max_abs}")
    print(f"  sum abs difference:           {sum_abs:,}")
    if nonzero_diff:
        print("  sample differing cells:")
        with pl.Config(tbl_rows=10):
            print(diff.filter(pl.col("absdiff") != 0)
                      .sort("absdiff", descending=True)
                      .head(10))
    print(f"\n  PASS (>99.9% identical): {frac > 0.999}")


def main() -> int:
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    ref, test = Path(sys.argv[1]), Path(sys.argv[2])
    print(f"REF  (H200 full):   {ref}")
    print(f"TEST (sharded RTX): {test}")
    compare_hits(ref, test)
    compare_counts(ref, test)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
