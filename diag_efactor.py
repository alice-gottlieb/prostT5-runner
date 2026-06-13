#!/usr/bin/env python3
"""Characterize the e-value ratio between full-DB and sharded runs on shared
(query,target) hits, and confirm bit scores match. The sharded run's evalue was
already multiplied by the shard count (6); this measures the TRUE factor.
Usage: diag_efactor.py REF_DIR TEST_DIR"""
import sys
from pathlib import Path
import polars as pl

ref, test = Path(sys.argv[1]), Path(sys.argv[2])


def load(d):
    return (pl.scan_csv(d / "hits.tsv", separator="\t")
            .select(["query", "target", "evalue", "bits"]).collect()
            .unique(["query", "target"]))


r = load(ref).rename({"evalue": "e_ref", "bits": "b_ref"})
t = load(test).rename({"evalue": "e_test", "bits": "b_test"})
j = r.join(t, on=["query", "target"], how="inner")
print(f"shared (query,target) hits: {j.height:,}")

bits_match = j.filter(pl.col("b_ref") == pl.col("b_test")).height
print(f"bit scores identical:       {bits_match:,} ({100*bits_match/j.height:.2f}%)")

# test evalue was multiplied by 6 already; true factor = (e_ref / (e_test/6)).
j = j.filter((pl.col("e_test") > 0) & (pl.col("e_ref") > 0)).with_columns(
    (pl.col("e_ref") / (pl.col("e_test") / 6)).alias("true_factor")
)
f = j["true_factor"]
print(f"\nTRUE e-value factor (e_ref / e_shard_raw):")
print(f"  median: {f.median():.4f}")
print(f"  mean:   {f.mean():.4f}")
print(f"  p05..p95: {f.quantile(0.05):.4f} .. {f.quantile(0.95):.4f}")
print(f"  (my code used 6; ratio my/true = {6/f.median():.3f})")
