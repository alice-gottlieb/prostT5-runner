#!/usr/bin/env python3
"""Among queries NOT hit by the max-seqs cap (<1000 hits in both runs), how well
do RTX-sharded and full-DB hits agree? Isolates the cap effect from the
sharding+e-rescale machinery. Usage: diag_noncapped.py REF_DIR TEST_DIR"""
import sys
from pathlib import Path
import polars as pl

ref, test = Path(sys.argv[1]), Path(sys.argv[2])
CAP = 1000


def pairs(d):
    return (pl.scan_csv(d / "hits.tsv", separator="\t")
            .select(["query", "target"]).unique().collect())


r, t = pairs(ref), pairs(test)
rn = r.group_by("query").len().rename({"len": "rn"})
tn = t.group_by("query").len().rename({"len": "tn"})

# Queries present in both and below the cap in both.
both = rn.join(tn, on="query", how="inner")
noncapped = both.filter((pl.col("rn") < CAP) & (pl.col("tn") < CAP))
print(f"queries in both runs:            {both.height:,}")
print(f"non-capped in both (<{CAP} hits): {noncapped.height:,}")

keep = noncapped.select("query")
rp = r.join(keep, on="query", how="inner")
tp = t.join(keep, on="query", how="inner")
shared = rp.join(tp, on=["query", "target"], how="inner").height
ru, tu = rp.height, tp.height
union = ru + tu - shared
print(f"\nNon-capped (query,target) pairs:")
print(f"  ref pairs:   {ru:,}")
print(f"  test pairs:  {tu:,}")
print(f"  shared:      {shared:,}")
print(f"  ref-only:    {ru - shared:,}")
print(f"  test-only:   {tu - shared:,}")
print(f"  Jaccard:     {shared/union if union else 1.0:.6f}")

# Also report agreement on per-query hit counts for non-capped queries.
exact = noncapped.filter(pl.col("rn") == pl.col("tn")).height
print(f"\nNon-capped queries with identical hit COUNT: {exact:,} "
      f"({100*exact/noncapped.height:.2f}%)")
