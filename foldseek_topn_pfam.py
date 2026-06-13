#!/usr/bin/env python3
"""Foldseek-search 3Di codes from a top-N Pfam-hits TSV and aggregate by genome.

Workflow:
  1. Pull 3Di codes from a TSV like all_pfam_top50_hits.tsv (one query per row).
  2. Build a Foldseek query DB and run a GPU-accelerated `foldseek search`
     against a Foldseek 3Di target DB.
  3. Record every hit, parse the target's genome, group by genome, and count
     how many times each 3Di code is found in each genome.

Outputs:
  <out>/hits.tsv               Long-form TSV: one foldseek hit per row plus
                               query metadata and target genome.
  <out>/genome_counts.tsv      Wide count matrix: rows = genome, cols = query_id,
                               values = hit count.

Examples:
  # Small local CPU smoke test.
  uv run python foldseek_topn_pfam.py pf00198_top5.tsv \
      results_fs_only_full_fs_compare/foldseek_db/sequenceDB \
      foldseek_topn_pfam_pf00198_test_out \
      --foldseek ~/foldseek/bin/foldseek --gpu 0 --threads 2 --max-seqs 5

  # Test subset from an all-Pfam top-N table.
  uv run python foldseek_topn_pfam.py all_pfam_top10_hits.tsv \
      results_fs_only_full_fs_compare/foldseek_db/sequenceDB \
      foldseek_topn_pfam_test_out \
      --foldseek ~/foldseek/bin/foldseek --gpu 0 --threads 2 --max-seqs 5 --test

  # Server-style run with regular flushed progress messages.
  uv run python foldseek_topn_pfam.py all_pfam_top50_hits.tsv /path/to/targetDB \
      foldseek_topn_pfam_out \
      --foldseek ~/foldseek/bin/foldseek --threads 16 --gpu 1 \
      --split-memory-limit 40G --progress-interval 1000

Test subset (--test) keeps only the top 5 rows of the first 3 Pfam accessions.
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from build_merged_prefilter import build_merged_prefilter

SEQUENCE_DBTYPE = b"\x00\x00\x00\x00"
HEADER_DBTYPE = b"\x0c\x00\x00\x00"

GENOME_RE = re.compile(r"(GC[FA]_\d+\.\d+)")


def log(message: str) -> None:
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}", flush=True)


def fresh_dir(path: Path) -> Path:
    """Remove `path` if present, then create it empty. Returns the path."""
    shutil.rmtree(path, ignore_errors=True)
    path.mkdir(parents=True)
    return path


def run_foldseek(foldseek: str, subcommand: str, *cmd_args: object) -> None:
    """Run a foldseek subcommand, stringifying each argument and logging it."""
    cmd = [foldseek, subcommand, *(str(arg) for arg in cmd_args)]
    log(f"Running foldseek {subcommand}: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


@dataclass(frozen=True)
class Query:
    query_id: str
    pfam_accession: str
    query_genome: str
    target_name: str
    rank: str
    source_file: str
    three_di: str


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input_tsv", type=Path,
                   help="TSV with target_name, rank, 3di_sequence, source_file, and query_pfam_accession or domain_accession columns.")
    p.add_argument("target_db", type=Path,
                   help="Full Foldseek target DB prefix (e.g. /path/to/merged3diDB). "
                        "With --target-shards this is the DB the merged prefilter is "
                        "aligned against (so E-values are full-DB values).")
    p.add_argument("output_dir", type=Path,
                   help="Directory to write hits.tsv, genome_counts.tsv, query_db/, etc.")
    p.add_argument("--foldseek", default=str(Path.home() / "foldseek" / "bin" / "foldseek"),
                   help="Path to foldseek binary (default: ~/foldseek/bin/foldseek).")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--gpu", type=int, default=1, choices=(0, 1),
                   help="Pass --gpu to foldseek search (1 = GPU, default).")
    p.add_argument("--split-memory-limit", default="40G")
    p.add_argument("--split", type=int, default=None,
                   help="Split the target DB into N chunks (foldseek --split). "
                        "Use with --split-mode 0 to fit a large target DB into "
                        "limited GPU VRAM; foldseek corrects e-values to the full "
                        "DB size, so results match an unsplit search.")
    p.add_argument("--split-mode", type=int, default=None, choices=(0, 1, 2),
                   help="foldseek --split-mode (0: split target DB, 1: split "
                        "query DB, 2: auto). Pair --split-mode 0 with --split.")
    p.add_argument("--max-seqs", type=int, default=None,
                   help="Optional cap on hits per query (passed as foldseek --max-seqs).")
    p.add_argument("--target-shards", type=Path, default=None,
                   help="Directory of padded shard DBs (prefixes named *_pad) that "
                        "together partition the full target_db. Enables GPU search "
                        "on a small card: the GPU ungapped prefilter runs per shard "
                        "(each fits VRAM), and the merged global top hits are aligned "
                        "against the full target_db, so E-values are full-DB values.")
    p.add_argument("--merge-max-seqs", type=int, default=1000,
                   help="Keep at most this many prefilter hits per query when "
                        "merging shards (mirrors foldseek's prefilter --max-seqs).")
    p.add_argument("--target-genome-map", type=Path, default=None,
                   help="Optional TSV mapping target_id -> genome. First two columns used, header skipped.")
    p.add_argument("--test", action="store_true",
                   help="Only use top 5 rows of the first 3 Pfam accessions in the TSV.")
    p.add_argument("--keep-tmp", action="store_true",
                   help="Keep foldseek tmp search directory.")
    p.add_argument("--progress-interval", type=int, default=1000,
                   help="Flush progress every N input/output rows. Use 0 to disable periodic progress logs.")
    return p.parse_args()


def load_queries(tsv: Path, test_mode: bool, progress_interval: int) -> list[Query]:
    """Read the input TSV and turn each non-empty 3Di row into a Query."""
    queries: list[Query] = []
    test_pfams: list[str] = []
    per_pfam_kept: Counter = Counter()
    rows_seen = 0
    rows_with_3di = 0

    log(f"Reading input TSV: {tsv}")
    with tsv.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"target_name", "rank", "3di_sequence", "source_file"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            sys.exit(f"ERROR: input TSV is missing required columns: {sorted(missing)}")
        if not {"query_pfam_accession", "domain_accession"} & set(reader.fieldnames or []):
            sys.exit(
                "ERROR: input TSV is missing a Pfam accession column: "
                "expected query_pfam_accession or domain_accession."
            )

        for row in reader:
            rows_seen += 1
            seq = (row.get("3di_sequence") or "").strip()
            if not seq or seq == "NA":
                if progress_interval and rows_seen % progress_interval == 0:
                    log(
                        f"Scanned {rows_seen:,} input rows; "
                        f"kept {len(queries):,} queries with 3Di sequences."
                    )
                continue
            rows_with_3di += 1

            pfam = row.get("query_pfam_accession") or row.get("domain_accession") or "NA"

            if test_mode:
                if pfam not in test_pfams:
                    if len(test_pfams) >= 3:
                        continue
                    test_pfams.append(pfam)
                if per_pfam_kept[pfam] >= 5:
                    continue
                per_pfam_kept[pfam] += 1

            target = row["target_name"]
            rank = row["rank"]
            source_file = row["source_file"]
            query_genome = parse_genome(source_file)
            qid = re.sub(
                r"\s+",
                "_",
                f"{pfam}|{query_genome}|{target}|rank{rank}",
            )
            queries.append(Query(
                query_id=qid,
                pfam_accession=pfam,
                query_genome=query_genome,
                target_name=target,
                rank=rank,
                source_file=source_file,
                three_di=seq,
            ))

            if progress_interval and rows_seen % progress_interval == 0:
                log(
                    f"Scanned {rows_seen:,} input rows; "
                    f"found {rows_with_3di:,} rows with 3Di; "
                    f"kept {len(queries):,} queries."
                )

    if not queries:
        sys.exit("ERROR: no non-empty 3Di sequences found in input TSV.")
    log(
        f"Finished reading {rows_seen:,} input rows; "
        f"found {rows_with_3di:,} rows with 3Di; "
        f"kept {len(queries):,} queries."
    )
    return queries


def write_query_db(queries: list[Query], output_dir: Path,
                   progress_interval: int) -> tuple[Path, Path]:
    """Write a Foldseek-compatible query DB by hand and a metadata TSV.

    Mirrors the byte layout used by build_extracted_3di_query_db.py.
    """
    query_db_dir = fresh_dir(output_dir / "query_db")
    query_db = query_db_dir / "queryDB"
    metadata_tsv = output_dir / "query_metadata.tsv"

    seq_offset = 0
    header_offset = 0

    log(f"Writing Foldseek query DB for {len(queries):,} queries: {query_db}")
    with (
        query_db.open("wb") as db,
        Path(f"{query_db}_ss").open("wb") as ss_db,
        Path(f"{query_db}_h").open("wb") as header_db,
        Path(f"{query_db}_ss_h").open("wb") as ss_header_db,
        Path(f"{query_db}.index").open("w") as index,
        Path(f"{query_db}_ss.index").open("w") as ss_index,
        Path(f"{query_db}_h.index").open("w") as header_index,
        Path(f"{query_db}_ss_h.index").open("w") as ss_header_index,
        (query_db_dir / "queryDB.lookup").open("w") as lookup,
        metadata_tsv.open("w", newline="") as meta,
    ):
        writer = csv.writer(meta, delimiter="\t", lineterminator="\n")
        writer.writerow(["query_id", "pfam_accession", "query_genome",
                         "target_name", "rank", "source_file",
                         "three_di_sequence"])

        for i, q in enumerate(queries):
            seq_bytes = q.three_di.encode()
            header_bytes = f"{q.query_id}\n\0".encode()

            db.write(seq_bytes)
            ss_db.write(seq_bytes)
            header_db.write(header_bytes)
            ss_header_db.write(header_bytes)

            index.write(f"{i}\t{seq_offset}\t{len(seq_bytes)}\n")
            ss_index.write(f"{i}\t{seq_offset}\t{len(seq_bytes)}\n")
            header_index.write(f"{i}\t{header_offset}\t{len(header_bytes)}\n")
            ss_header_index.write(f"{i}\t{header_offset}\t{len(header_bytes)}\n")
            lookup.write(f"{i}\t{q.query_id}\t0\n")

            writer.writerow([q.query_id, q.pfam_accession, q.query_genome,
                             q.target_name, q.rank, q.source_file,
                             q.three_di])

            seq_offset += len(seq_bytes)
            header_offset += len(header_bytes)

            if progress_interval and (i + 1) % progress_interval == 0:
                log(f"Wrote {i + 1:,}/{len(queries):,} query DB entries.")

    Path(f"{query_db}.dbtype").write_bytes(SEQUENCE_DBTYPE)
    Path(f"{query_db}_ss.dbtype").write_bytes(SEQUENCE_DBTYPE)
    Path(f"{query_db}_h.dbtype").write_bytes(HEADER_DBTYPE)
    Path(f"{query_db}_ss_h.dbtype").write_bytes(HEADER_DBTYPE)

    log(f"Finished writing query DB and metadata TSV: {metadata_tsv}")
    return query_db, metadata_tsv


FOLDSEEK_COLUMNS = ["query", "target", "fident", "alnlen", "mismatch",
                    "gapopen", "qstart", "qend", "tstart", "tend",
                    "evalue", "bits"]


def run_foldseek_search(args: argparse.Namespace, query_db: Path,
                        output_dir: Path, target_db: Path,
                        tag: str = "") -> Path:
    result_prefix = output_dir / f"result{tag}"
    tmp_dir = fresh_dir(output_dir / f"tmp_search{tag}")

    options = [
        "--threads", args.threads,
        "--gpu", args.gpu,
        "--split-memory-limit", args.split_memory_limit,
    ]
    if args.split is not None:
        options += ["--split", args.split]
    if args.split_mode is not None:
        options += ["--split-mode", args.split_mode]
    if args.max_seqs is not None:
        options += ["--max-seqs", args.max_seqs]

    run_foldseek(args.foldseek, "search",
                 query_db, target_db, result_prefix, tmp_dir, *options)

    if not Path(f"{result_prefix}.index").exists():
        sys.exit(f"ERROR: foldseek search produced no {result_prefix}.index")
    log(f"Foldseek search finished: {result_prefix}")
    return result_prefix


def run_convertalis(foldseek: str, query_db: Path, target_db: Path,
                    result_prefix: Path, out_tsv: Path, threads: int,
                    extra_format: bool) -> None:
    options = ["--threads", threads]
    if extra_format:
        options += ["--format-output", ",".join(FOLDSEEK_COLUMNS)]
    run_foldseek(foldseek, "convertalis",
                 query_db, target_db, result_prefix, out_tsv, *options)
    log(f"Foldseek convertalis finished: {out_tsv}")


# Parameters copied verbatim from foldseek search's own GPU pipeline so the
# sharded run reproduces it exactly: the ungapped prefilter call (captured from a
# `--gpu 1` search) and the structurealign call (from foldseek's structuresearch.sh).
UNGAPPED_PREFILTER_PAR = [
    "--sub-mat", "aa:3di.out,nucl:3di.out", "-c", 0, "-e", "1.79769e+308",
    "--cov-mode", 0, "--comp-bias-corr", 1, "--comp-bias-corr-scale", 0.15,
    "--min-ungapped-score", 30, "--max-seqs", 1000, "--db-load-mode", 0,
    "--gpu", 1, "-v", 3,
]
STRUCTUREALIGN_PAR = [
    "--tmscore-threshold", 0, "--tmscore-threshold-mode", 0, "--lddt-threshold", 0,
    "--sort-by-structure-bits", 1, "--alignment-type", 2, "--exact-tmscore", 0,
    "--sub-mat", "aa:3di.out,nucl:3di.out", "-a", 0, "--alignment-mode", 3,
    "--alignment-output-mode", 0, "--wrapped-scoring", 0, "-e", 10,
    "--min-seq-id", 0, "--min-aln-len", 0, "--seq-id-mode", 0, "--alt-ali", 0,
    "-c", 0, "--cov-mode", 0, "--max-seq-len", 65535, "--comp-bias-corr", 1,
    "--comp-bias-corr-scale", 0.5, "--max-rejected", 2147483647,
    "--max-accept", 2147483647, "--add-self-matches", 0, "--db-load-mode", 0,
    "--score-bias", 0, "--realign", 0, "--gap-open", "aa:10,nucl:10",
    "--gap-extend", "aa:1,nucl:1", "--zdrop", 40, "-v", 3,
]


def discover_shards(shards_dir: Path) -> list[Path]:
    """Padded shard DB prefixes in `shards_dir` (files named <prefix>_pad)."""
    prefixes = sorted(p.with_suffix("") for p in shards_dir.glob("*_pad.dbtype"))
    if not prefixes:
        sys.exit(f"ERROR: no padded shard DBs (*_pad) found in {shards_dir}")
    return prefixes


def run_sharded_prefilter_merge(args: argparse.Namespace, query_db: Path,
                                output_dir: Path) -> Path:
    """Reproduce a full-DB GPU search by sharding only the prefilter.

    foldseek's GPU search loads the whole padded target DB into VRAM, which a
    small card cannot hold. Instead we run the GPU ungapped prefilter against
    each padded shard (each fits VRAM), resolve target keys to accessions per
    shard (createtsv), then merge into one global top-N prefilter in the full-DB
    keyspace. We align that with structurealign against the WHOLE DB, so the
    E-values are the true full-DB values, and convertalis. This mirrors
    foldseek's own ungappedprefilter -> structurealign -> convert pipeline.
    """
    shards = discover_shards(args.target_shards)
    full_db = args.target_db
    threads = args.threads
    prefilter_dir = fresh_dir(output_dir / "prefilter")
    log(f"Sharded prefilter+merge over {len(shards)} shards against {full_db}.")

    tsvs: list[str] = []
    for shard in shards:
        pref = prefilter_dir / f"pref_{shard.name}"
        log(f"GPU ungapped prefilter on shard {shard.name}")
        run_foldseek(args.foldseek, "ungappedprefilter",
                     f"{query_db}_ss", f"{shard}_ss", pref,
                     *UNGAPPED_PREFILTER_PAR, "--threads", threads)
        tsv = prefilter_dir / f"pref_{shard.name}.tsv"
        run_foldseek(args.foldseek, "createtsv", query_db, shard, pref, tsv,
                     "--threads", threads)
        tsvs.append(str(tsv))

    pref_merged = prefilter_dir / "pref_merged"
    log("Merging shard prefilters into the full-DB keyspace.")
    build_merged_prefilter(str(query_db), str(full_db), str(pref_merged),
                           args.merge_max_seqs, tsvs, log=log)

    aln = output_dir / "aln"
    log("structurealign against the full DB (true full-DB E-values).")
    run_foldseek(args.foldseek, "structurealign", query_db, full_db,
                 pref_merged, aln, *STRUCTUREALIGN_PAR, "--threads", threads)

    results_full = output_dir / "results.full.tsv"
    run_convertalis(args.foldseek, query_db, full_db, aln, results_full,
                    threads, extra_format=True)

    if not args.keep_tmp:
        shutil.rmtree(prefilter_dir, ignore_errors=True)
    return results_full


def load_target_genome_map(path: Path | None) -> dict[str, list[str]]:
    """Load a target_id -> [genome, ...] map (one row per target/genome pair).

    A single target (e.g. an NCBI WP_ multispecies protein) can live in many
    genomes, so each target maps to a list of genomes rather than just one. A
    hit on that target then counts toward every genome it appears in.
    """
    if path is None:
        log("No target genome map provided; parsing genomes from target IDs.")
        return {}
    mapping: dict[str, list[str]] = defaultdict(list)
    log(f"Loading target genome map: {path}")
    with path.open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader, None)  # header
        for row in reader:
            if len(row) >= 2 and row[0]:
                mapping[row[0]].append(row[1])
    log(f"Loaded genome mappings for {len(mapping):,} targets.")
    return dict(mapping)


def parse_genome(text: str) -> str:
    """Pull a single genome accession (GCF_/GCA_) out of a string, or 'NA'."""
    m = GENOME_RE.search(text)
    return m.group(1) if m else "NA"


def genomes_for_target(target: str, mapping: dict[str, list[str]]) -> list[str]:
    """Every genome a target belongs to: from the map, else parsed from the id."""
    if target in mapping:
        return mapping[target]
    return [parse_genome(target)]


def write_hits_and_counts(results_tsv: Path, queries: list[Query],
                          target_map: dict[str, str], output_dir: Path,
                          progress_interval: int) -> None:
    hits_tsv = output_dir / "hits.tsv"
    counts_tsv = output_dir / "genome_counts.tsv"

    query_meta = {q.query_id: q for q in queries}

    # genome -> query_id -> count
    counts: dict[str, Counter] = defaultdict(Counter)
    all_query_ids: list[str] = [q.query_id for q in queries]
    hits_seen = 0

    log(f"Aggregating Foldseek results: {results_tsv}")
    with results_tsv.open() as src, hits_tsv.open("w", newline="") as dst:
        writer = csv.writer(dst, delimiter="\t", lineterminator="\n")
        writer.writerow(FOLDSEEK_COLUMNS + ["target_genome", "pfam_accession",
                                            "query_genome",
                                            "query_target_name", "query_rank",
                                            "query_source_file"])
        for line in src:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(FOLDSEEK_COLUMNS):
                continue
            query_id, target = parts[0], parts[1]
            genomes = genomes_for_target(target, target_map)
            q = query_meta.get(query_id)
            writer.writerow(parts + [
                ";".join(genomes),
                q.pfam_accession if q else "NA",
                q.query_genome if q else "NA",
                q.target_name if q else "NA",
                q.rank if q else "NA",
                q.source_file if q else "NA",
            ])
            # A multispecies target lives in many genomes; count the hit in each.
            for genome in genomes:
                counts[genome][query_id] += 1
            hits_seen += 1

            if progress_interval and hits_seen % progress_interval == 0:
                log(
                    f"Aggregated {hits_seen:,} hits across "
                    f"{len(counts):,} target genomes."
                )

    genomes = sorted(counts)
    log(f"Writing genome count matrix for {len(genomes):,} genomes.")
    with counts_tsv.open("w", newline="") as dst:
        writer = csv.writer(dst, delimiter="\t", lineterminator="\n")
        writer.writerow(["genome"] + all_query_ids)
        for genome in genomes:
            row = [genome] + [counts[genome].get(qid, 0) for qid in all_query_ids]
            writer.writerow(row)

    log(f"Wrote {hits_tsv} ({hits_seen:,} hits).")
    log(f"Wrote {counts_tsv}  ({len(genomes):,} genomes x {len(all_query_ids):,} queries).")


def main() -> int:
    try:
        sys.stdout.reconfigure(line_buffering=True)
        sys.stderr.reconfigure(line_buffering=True)
    except AttributeError:
        pass

    args = parse_args()
    log("Starting foldseek_topn_pfam.py")
    log(f"Arguments: input_tsv={args.input_tsv}, target_db={args.target_db}, output_dir={args.output_dir}")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    queries = load_queries(args.input_tsv, args.test, args.progress_interval)
    log(f"Loaded {len(queries):,} 3Di queries"
        f" ({'TEST subset' if args.test else 'full input'}).")

    query_db, _ = write_query_db(queries, args.output_dir, args.progress_interval)
    log(f"Built query DB: {query_db}")

    if args.target_shards:
        results_full = run_sharded_prefilter_merge(args, query_db, args.output_dir)
    else:
        result_prefix = run_foldseek_search(args, query_db, args.output_dir,
                                            args.target_db)

        # Two convertalis flavours run in parallel: default columns and a richer
        # format we actually parse. Both use foldseek's own thread pool, so we keep
        # this to two workers to avoid oversubscribing the GPU/CPU.
        results_default = args.output_dir / "results.tsv"
        results_full = args.output_dir / "results.full.tsv"

        log("Starting Foldseek convertalis jobs.")
        with ThreadPoolExecutor(max_workers=2) as pool:
            futures = [
                pool.submit(run_convertalis, args.foldseek, query_db, args.target_db,
                            result_prefix, results_default, args.threads, False),
                pool.submit(run_convertalis, args.foldseek, query_db, args.target_db,
                            result_prefix, results_full, args.threads, True),
            ]
            for f in futures:
                f.result()
        log("Finished Foldseek convertalis jobs.")

    target_map = load_target_genome_map(args.target_genome_map)
    write_hits_and_counts(results_full, queries, target_map, args.output_dir,
                          args.progress_interval)

    if not args.keep_tmp:
        shutil.rmtree(args.output_dir / "tmp_search", ignore_errors=True)
        log(f"Removed temporary search directory: {args.output_dir / 'tmp_search'}")

    log("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
