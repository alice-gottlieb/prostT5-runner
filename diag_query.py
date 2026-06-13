#!/usr/bin/env python3
"""Isolate truly non-capped queries (raw hit ROWS < 1000 in both) and inspect one
query's hits side by side. Usage: diag_query.py REF_DIR TEST_DIR"""
import sys
from pathlib import Path
import polars as pl

ref, test = Path(sys.argv[1]), Path(sys.argv[2])
CAP = 1000


def rows(d):
    return pl.scan_csv(d / "hits.tsv", separator="\t").collect()


R, T = rows(ref), rows(test)
rc = R.group_by("query").len().rename({"len": "rrows"})
tc = T.group_by("query").len().rename({"len": "trows"})
both = rc.join(tc, on="query", how="inner")
noncap = both.filter((pl.col("rrows") < CAP) & (pl.col("trows") < CAP))
print(f"queries in both:               {both.height:,}")
print(f"truly non-capped (rows<{CAP}):  {noncap.height:,}")

# Jaccard on unique pairs restricted to truly non-capped queries.
keep = noncap.select("query")
rp = R.select(["query", "target"]).unique().join(keep, on="query", how="inner")
tp = T.select(["query", "target"]).unique().join(keep, on="query", how="inner")
shared = rp.join(tp, on=["query", "target"], how="inner").height
union = rp.height + tp.height - shared
print(f"non-capped Jaccard (pairs):    {shared/union if union else 1.0:.6f} "
      f"(ref {rp.height:,}, test {tp.height:,}, shared {shared:,})")

# Inspect one mid-sized non-capped query.
sample = noncap.filter((pl.col("rrows") > 100) & (pl.col("rrows") < 300))
if sample.height:
    q = sample["query"][0]
    print(f"\n=== sample non-capped query: {q}")
    for name, D in [("REF", R), ("TEST", T)]:
        h = (D.filter(pl.col("query") == q)
             .select(["target", "evalue", "bits"])
             .sort("bits", descending=True))
        print(f"  {name}: {h.height} rows, "
              f"{h.select(pl.col('target').n_unique()).item()} unique targets, "
              f"evalue max {h['evalue'].max():.3g}")
        with pl.Config(tbl_rows=6):
            print(h.head(6))
    rt = set(R.filter(pl.col("query") == q)["target"].to_list())
    tt = set(T.filter(pl.col("query") == q)["target"].to_list())
    print(f"  shared targets: {len(rt & tt)}, ref-only: {len(rt - tt)}, "
          f"test-only: {len(tt - rt)}")
