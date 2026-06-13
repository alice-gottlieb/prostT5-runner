#!/usr/bin/env python3
"""Merge per-shard ungappedprefilter results into one prefilter DB in the FULL-DB
keyspace, keeping the global top-N hits per query by ungapped score.

Each shard is GPU-prefiltered separately and re-keyed to a local keyspace, so we
resolve keys to accessions per shard (via `foldseek createtsv`), merge by
accession, then re-key targets to the full DB so `foldseek structurealign` can
align against the whole DB and compute true full-DB E-values.

Inputs are createtsv dumps of the per-shard prefilter results, each row:
    query_accession <tab> target_accession <tab> ungapped_score <tab> diagonal

Usage:
  build_merged_prefilter.py --query-db Q --full-db FULL --out PREF --max-seqs 1000 pref_*.tsv
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

PREFILTER_DBTYPE = (7).to_bytes(4, "little")


def load_lookup(prefix: str) -> dict[str, str]:
    """accession -> internal key, from a foldseek <prefix>.lookup file."""
    acc_to_key: dict[str, str] = {}
    with open(f"{prefix}.lookup") as handle:
        for line in handle:
            key, acc, _ = line.rstrip("\n").split("\t")
            acc_to_key[acc] = key
    return acc_to_key


def build_merged_prefilter(query_db: str, full_db: str, out: str,
                           max_seqs: int, tsvs: list[str],
                           log=print) -> None:
    """Merge per-shard createtsv prefilter dumps into one full-DB-keyspace
    prefilter DB, keeping the global top-`max_seqs` hits per query by score."""
    qacc_to_key = load_lookup(query_db)

    # Gather per query the (score, target_accession, diagonal) of every shard hit.
    per_query: dict[str, list[tuple[int, str, str]]] = defaultdict(list)
    needed_targets: set[str] = set()
    for tsv in tsvs:
        with open(tsv) as handle:
            for line in handle:
                qacc, tacc, score, diag = line.rstrip("\n").split("\t")
                per_query[qacc].append((int(score), tacc, diag))
                needed_targets.add(tacc)
    log(f"queries with hits: {len(per_query):,}; "
        f"unique targets: {len(needed_targets):,}")

    # accession -> full-DB key, only for targets we actually need.
    tacc_to_key: dict[str, str] = {}
    with open(f"{full_db}.lookup") as handle:
        for line in handle:
            key, acc, _ = line.rstrip("\n").split("\t")
            if acc in needed_targets:
                tacc_to_key[acc] = key
    log(f"resolved {len(tacc_to_key):,} target keys in full DB")

    # One entry per query (all queries, empty if no hits), keyed by query key,
    # keeping the global top-N by ungapped score.
    entries: list[tuple[int, bytes]] = []
    for qacc, qkey in qacc_to_key.items():
        hits = per_query.get(qacc, [])
        hits.sort(key=lambda h: h[0], reverse=True)
        lines = []
        for score, tacc, diag in hits[:max_seqs]:
            tkey = tacc_to_key.get(tacc)
            if tkey is not None:
                lines.append(f"{tkey}\t{score}\t{diag}\n")
        entries.append((int(qkey), "".join(lines).encode() + b"\x00"))

    entries.sort()
    out_path = Path(out)
    offset = 0
    with out_path.open("wb") as data, Path(f"{out}.index").open("w") as index:
        for qkey, entry in entries:
            data.write(entry)
            index.write(f"{qkey}\t{offset}\t{len(entry)}\n")
            offset += len(entry)
    Path(f"{out}.dbtype").write_bytes(PREFILTER_DBTYPE)
    log(f"wrote {len(entries):,} prefilter entries to {out}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--query-db", required=True, help="query DB prefix (.lookup)")
    ap.add_argument("--full-db", required=True, help="full target DB prefix (.lookup)")
    ap.add_argument("--out", required=True, help="output prefilter DB prefix")
    ap.add_argument("--max-seqs", type=int, default=1000)
    ap.add_argument("tsvs", nargs="+", help="per-shard createtsv prefilter dumps")
    args = ap.parse_args()
    build_merged_prefilter(args.query_db, args.full_db, args.out,
                           args.max_seqs, args.tsvs)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
