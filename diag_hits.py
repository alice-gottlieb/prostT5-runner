#!/usr/bin/env python3
"""Diagnose RTX-sharded vs full-DB hit divergence: how cap-bound are queries,
and how does one divergent query differ?  Usage: diag_hits.py REF_DIR TEST_DIR"""
import sys
from pathlib import Path
import polars as pl

ref, test = Path(sys.argv[1]), Path(sys.argv[2])
probe = sys.argv[3] if len(sys.argv) > 3 else None

for name, d in [("REF (full)", ref), ("TEST (sharded)", test)]:
    n = (pl.scan_csv(d / "hits.tsv", separator="\t")
         .group_by("query").agg(pl.len().alias("n")).collect())["n"]
    print(f"\n{name}: {d}")
    print(f"  queries:            {n.len():,}")
    print(f"  total hits:         {n.sum():,}")
    print(f"  mean hits/query:    {n.mean():.1f}")
    print(f"  median hits/query:  {n.median():.0f}")
    print(f"  max hits/query:     {n.max():,}")
    for thr in (1000, 999, 500, 100):
        print(f"  queries with >= {thr:>4} hits: {n.filter(n >= thr).len():,} "
              f"({100*n.filter(n >= thr).len()/n.len():.1f}%)")

if probe:
    print(f"\n=== probe query: {probe} ===")
    for name, d in [("REF", ref), ("TEST", test)]:
        df = (pl.scan_csv(d / "hits.tsv", separator="\t")
              .filter(pl.col("query") == probe)
              .select(["target", "evalue", "bits"]).collect())
        print(f"  {name}: {df.height} hits; "
              f"bits range [{df['bits'].min()}, {df['bits'].max()}]" if df.height
              else f"  {name}: 0 hits")
