#!/usr/bin/env python3
"""Output-integrity tests for foldseek_topn_pfam.py.

These prove that every cell in hits.tsv / genome_counts.tsv traces back to the
correct input row: the right pfam, the right genome, no shifted or mis-copied
columns. Cases are built deliberately to catch mis-joins (wrong query_id ->
metadata), positional column copy bugs, and byte-offset accumulation bugs in the
hand-written Foldseek query DB.

The data-logic groups (A-D) need no foldseek binary and run anywhere. Group E is
an end-to-end test that shells out to the real foldseek binary; it is skipped
when the binary or target DB is absent.

Run:
    uv run python test_foldseek_topn_pfam.py
    FOLDSEEK=~/foldseek/bin/foldseek uv run python test_foldseek_topn_pfam.py

Exits 0 and prints ALL TESTS PASSED only when every non-skipped check passes;
exits 1 on the first failed check.
"""

from __future__ import annotations

import csv
import os
import subprocess
import sys
import tempfile
from pathlib import Path

from foldseek_topn_pfam import (
    FOLDSEEK_COLUMNS,
    HEADER_DBTYPE,
    SEQUENCE_DBTYPE,
    Query,
    load_queries,
    load_target_genome_map,
    target_to_genome,
    write_hits_and_counts,
    write_query_db,
)

REPO = Path(__file__).resolve().parent


# --------------------------------------------------------------------------- #
# Tiny harness, matching test_merge.py's style.
# --------------------------------------------------------------------------- #
def section(title: str) -> None:
    print(f"\n=== {title} ===")


def fail(msg: str) -> None:
    print(f"FAIL: {msg}", file=sys.stderr)
    sys.exit(1)


def check(cond: bool, msg: str) -> None:
    if not cond:
        fail(msg)
    print(f"  ok: {msg}")


def write_tsv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def read_tsv(path: Path) -> tuple[list[str], list[list[str]]]:
    with path.open(newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    return rows[0], rows[1:]


# --------------------------------------------------------------------------- #
# Group A: load_queries (input row -> Query)
# --------------------------------------------------------------------------- #
# Column order for the domain_accession variant (matches pf00198_top5.tsv).
DOMAIN_HEADER = [
    "rank", "target_name", "domain_name", "domain_accession",
    "3di_sequence", "source_file",
]


def domain_row(rank, target, pfam, three_di, source) -> list[str]:
    return [rank, target, "name", pfam, three_di, source]


def test_a_load_queries(tmp: Path) -> None:
    section("A: load_queries")

    # A1 pfam precedence: query_pfam_accession wins over domain_accession.
    both = tmp / "both.tsv"
    write_tsv(
        both,
        ["query_pfam_accession", "domain_accession", "rank", "target_name",
         "3di_sequence", "source_file"],
        [["PF99999.1", "PF00001.1", "1", "T1", "ABCD", "GCF_000000001.1.x"]],
    )
    q = load_queries(both, test_mode=False, progress_interval=0)
    check(len(q) == 1 and q[0].pfam_accession == "PF99999.1",
          "A1 query_pfam_accession takes precedence over domain_accession")

    # A1b only domain_accession present -> used as fallback.
    only_dom = tmp / "only_dom.tsv"
    write_tsv(only_dom, DOMAIN_HEADER,
              [domain_row("1", "T1", "PF00198.29", "ABCD", "GCF_000005845.2.x")])
    q = load_queries(only_dom, test_mode=False, progress_interval=0)
    check(q[0].pfam_accession == "PF00198.29",
          "A1b falls back to domain_accession when no query_pfam_accession")

    # A1c both pfam columns present but empty on the row -> "NA".
    empty_pfam = tmp / "empty_pfam.tsv"
    write_tsv(
        empty_pfam,
        ["query_pfam_accession", "domain_accession", "rank", "target_name",
         "3di_sequence", "source_file"],
        [["", "", "1", "T1", "ABCD", "GCF_000000001.1.x"]],
    )
    q = load_queries(empty_pfam, test_mode=False, progress_interval=0)
    check(q[0].pfam_accession == "NA",
          "A1c missing both pfam values yields NA")

    # A2 skip empty / whitespace / NA 3Di.
    skip = tmp / "skip.tsv"
    write_tsv(skip, DOMAIN_HEADER, [
        domain_row("1", "Keep", "PF1", "ABCD", "GCF_000000001.1.x"),
        domain_row("2", "Empty", "PF1", "", "GCF_000000001.1.x"),
        domain_row("3", "Space", "PF1", "   ", "GCF_000000001.1.x"),
        domain_row("4", "Na", "PF1", "NA", "GCF_000000001.1.x"),
    ])
    q = load_queries(skip, test_mode=False, progress_interval=0)
    check([x.target_name for x in q] == ["Keep"],
          "A2 empty / whitespace / NA 3Di rows are dropped")

    # A3 query_id format + whitespace collapse.
    ws = tmp / "ws.tsv"
    write_tsv(ws, DOMAIN_HEADER,
              [domain_row("1", "T 1", "PF1", "ABCD", "GCF_000005845.2.x")])
    q = load_queries(ws, test_mode=False, progress_interval=0)
    check(q[0].query_id == "PF1|GCF_000005845.2|T_1|rank1",
          "A3 query_id format with whitespace collapsed to underscore")
    check(" " not in q[0].query_id, "A3 no whitespace survives in query_id")

    # A4 query_genome parsed from source_file.
    gen = tmp / "gen.tsv"
    write_tsv(gen, DOMAIN_HEADER, [
        domain_row("1", "T1", "PF1", "ABCD", "GCF_000005845.2.with_3di.domtblout"),
        domain_row("2", "T2", "PF1", "ABCD", "no_accession_here.txt"),
    ])
    q = load_queries(gen, test_mode=False, progress_interval=0)
    check(q[0].query_genome == "GCF_000005845.2",
          "A4 query_genome parsed from source_file accession")
    check(q[1].query_genome == "NA",
          "A4 source_file without accession yields NA query_genome")

    # A5 --test subset: first 3 distinct pfams, <=5 rows each.
    sub = tmp / "sub.tsv"
    rows = []
    for pfam in ("PFa", "PFb", "PFc", "PFd"):
        for r in range(1, 7):  # 6 rows each
            rows.append(domain_row(str(r), f"T_{pfam}_{r}", pfam, "ABCD",
                                   "GCF_000000001.1.x"))
    write_tsv(sub, DOMAIN_HEADER, rows)
    q = load_queries(sub, test_mode=True, progress_interval=0)
    kept = {}
    for x in q:
        kept[x.pfam_accession] = kept.get(x.pfam_accession, 0) + 1
    check(len(q) == 15, "A5 --test keeps exactly 15 queries (3 pfams x 5)")
    check(kept == {"PFa": 5, "PFb": 5, "PFc": 5},
          "A5 --test keeps top 5 of the first 3 pfams only")

    # A6 missing-column errors -> SystemExit.
    miss_req = tmp / "miss_req.tsv"
    write_tsv(miss_req, ["rank", "target_name", "domain_accession",
                         "3di_sequence"],  # no source_file
              [["1", "T1", "PF1", "ABCD"]])
    try:
        load_queries(miss_req, test_mode=False, progress_interval=0)
        fail("A6 expected SystemExit on missing required column")
    except SystemExit:
        print("  ok: A6 missing required column raises SystemExit")

    miss_pfam = tmp / "miss_pfam.tsv"
    write_tsv(miss_pfam, ["rank", "target_name", "3di_sequence", "source_file"],
              [["1", "T1", "ABCD", "GCF_000000001.1.x"]])
    try:
        load_queries(miss_pfam, test_mode=False, progress_interval=0)
        fail("A6 expected SystemExit on missing pfam column")
    except SystemExit:
        print("  ok: A6 missing both pfam columns raises SystemExit")


# --------------------------------------------------------------------------- #
# Group B: write_query_db byte-offset / metadata integrity (round-trip)
# --------------------------------------------------------------------------- #
def load_real_3di(n: int) -> list[str]:
    """Pull the first `n` non-empty, full 3Di strings from all_pfam_top10_hits.tsv."""
    src = REPO / "all_pfam_top10_hits.tsv"
    if not src.exists():
        fail(f"B requires real 3Di source file: {src}")
    seqs: list[str] = []
    with src.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            s = (row.get("3di_sequence") or "").strip()
            if s and s != "NA":
                seqs.append(s)
                if len(seqs) >= n:
                    break
    if len(seqs) < n:
        fail(f"B needs {n} real 3Di strings, found only {len(seqs)}")
    return seqs


def parse_index(path: Path) -> list[tuple[int, int, int]]:
    out = []
    for line in path.read_text().splitlines():
        idx, offset, length = line.split("\t")
        out.append((int(idx), int(offset), int(length)))
    return out


def test_b_query_db(tmp: Path) -> None:
    section("B: write_query_db byte/metadata integrity (>=10 real 3Di)")

    real = load_real_3di(12)
    check(len(real) >= 10, "B uses >=10 real full 3Di strings")
    lengths = {len(s) for s in real}
    check(len(lengths) > 1, "B 3Di strings have varying lengths (exercises offsets)")

    queries = [
        Query(
            query_id=f"PF{i:05d}.1|GCF_000000{i:03d}.1|T{i}|rank{i}",
            pfam_accession=f"PF{i:05d}.1",
            query_genome=f"GCF_000000{i:03d}.1",
            target_name=f"T{i}",
            rank=str(i),
            source_file=f"GCF_000000{i:03d}.1.x",
            three_di=seq,
        )
        for i, seq in enumerate(real)
    ]

    out = tmp / "B"
    out.mkdir()
    query_db, metadata_tsv = write_query_db(queries, out, progress_interval=0)

    db_bytes = query_db.read_bytes()
    ss_bytes = Path(f"{query_db}_ss").read_bytes()
    h_bytes = Path(f"{query_db}_h").read_bytes()
    index = parse_index(Path(f"{query_db}.index"))
    ss_index = parse_index(Path(f"{query_db}_ss.index"))
    h_index = parse_index(Path(f"{query_db}_h.index"))

    # B1 index round-trip on the main sequence DB.
    ok = True
    for (idx, offset, length), q in zip(index, queries):
        sliced = db_bytes[offset:offset + length].decode()
        if idx != queries.index(q) or sliced != q.three_di:
            ok = False
            break
    check(len(index) == len(queries) and ok,
          "B1 every queryDB.index entry slices back to the right 3Di")

    # B2 _ss mirrors main exactly.
    check(ss_bytes == db_bytes and ss_index == index,
          "B2 queryDB_ss bytes and index mirror queryDB")

    # B3 header round-trip.
    ok = True
    for (idx, offset, length), q in zip(h_index, queries):
        sliced = h_bytes[offset:offset + length].decode()
        if sliced != f"{q.query_id}\n\0":
            ok = False
            break
    check(len(h_index) == len(queries) and ok,
          "B3 every header slices back to '{query_id}\\n\\0'")

    # B4 lookup file.
    lookup_lines = (out / "query_db" / "queryDB.lookup").read_text().splitlines()
    ok = all(
        lookup_lines[i] == f"{i}\t{q.query_id}\t0"
        for i, q in enumerate(queries)
    )
    check(len(lookup_lines) == len(queries) and ok,
          "B4 queryDB.lookup line i == 'i\\t{query_id}\\t0'")

    # B5 metadata alignment, field-for-field.
    header, data = read_tsv(metadata_tsv)
    check(header == ["query_id", "pfam_accession", "query_genome",
                     "target_name", "rank", "source_file", "three_di_sequence"],
          "B5 metadata header is in the expected order")
    ok = True
    for row, q in zip(data, queries):
        if row != [q.query_id, q.pfam_accession, q.query_genome, q.target_name,
                   q.rank, q.source_file, q.three_di]:
            ok = False
            break
    check(len(data) == len(queries) and ok,
          "B5 each metadata row matches its query field-for-field")

    # B6 dbtype bytes.
    check(Path(f"{query_db}.dbtype").read_bytes() == SEQUENCE_DBTYPE
          and Path(f"{query_db}_ss.dbtype").read_bytes() == SEQUENCE_DBTYPE,
          "B6 sequence dbtype bytes correct for main and _ss")
    check(Path(f"{query_db}_h.dbtype").read_bytes() == HEADER_DBTYPE
          and Path(f"{query_db}_ss_h.dbtype").read_bytes() == HEADER_DBTYPE,
          "B6 header dbtype bytes correct for _h and _ss_h")


# --------------------------------------------------------------------------- #
# Group C: target_to_genome / load_target_genome_map
# --------------------------------------------------------------------------- #
def test_c_genome_resolution(tmp: Path) -> None:
    section("C: target_to_genome / load_target_genome_map")

    # C1 map precedence: map wins even if target also contains a different accession.
    mapping = {"GCF_111111111.1_prot": "GCF_999999999.9"}
    check(target_to_genome("GCF_111111111.1_prot", mapping) == "GCF_999999999.9",
          "C1 map value wins over embedded accession in the target id")

    # C2 regex fallback.
    check(target_to_genome("WP_x_GCA_000123456.2_y", {}) == "GCA_000123456.2",
          "C2 accession parsed from target when not in map")
    check(target_to_genome("NP_459722.1", {}) == "NA",
          "C2 target without accession yields NA")

    # C3 map loading: header skipped, two cols, empty-first-col ignored.
    mp = tmp / "map.tsv"
    write_tsv(mp, ["target_id", "genome", "extra"], [
        ["t1", "GCF_000000001.1", "ignored"],
        ["", "GCF_000000002.2", "x"],          # empty first col -> skipped
        ["t3", "GCF_000000003.3", "y"],
    ])
    loaded = load_target_genome_map(mp)
    check(loaded == {"t1": "GCF_000000001.1", "t3": "GCF_000000003.3"},
          "C3 map loads two columns, skips header and empty-key rows")
    check(load_target_genome_map(None) == {},
          "C3 None path returns empty map")


# --------------------------------------------------------------------------- #
# Group D: write_hits_and_counts (the core join + exact output)
# --------------------------------------------------------------------------- #
def fs_line(query_id: str, target: str) -> str:
    """A synthetic 12-column foldseek convertalis line."""
    cols = [query_id, target, "0.5", "100", "10", "0",
            "1", "100", "1", "100", "1e-10", "200"]
    return "\t".join(cols)


def make_query(pfam, genome, target, rank) -> Query:
    qid = f"{pfam}|{genome}|{target}|rank{rank}"
    return Query(query_id=qid, pfam_accession=pfam, query_genome=genome,
                 target_name=target, rank=rank, source_file=f"{genome}.x",
                 three_di="ABCD")


def test_d_hits_and_counts(tmp: Path) -> None:
    section("D: write_hits_and_counts join + exact output")

    out = tmp / "D"
    out.mkdir()

    # Three queries with deliberately similar ids but distinct metadata, plus a
    # fourth that will get zero hits.
    q1 = make_query("PF00001.1", "GCF_000000001.1", "TGT_A", "1")
    q2 = make_query("PF00002.2", "GCF_000000002.2", "TGT_A", "2")  # same target diff pfam/genome
    q3 = make_query("PF00003.3", "GCF_000000003.3", "TGT_B", "1")
    q4 = make_query("PF00004.4", "GCF_000000004.4", "TGT_C", "1")  # zero hits
    queries = [q1, q2, q3, q4]

    # Targets resolve to genomes via a map so target_genome != query_genome,
    # proving the two genome columns are never crossed.
    target_map = {
        "hitX_GCF_000777777.7": "GCF_000777777.7",
        "hitY_GCF_000888888.8": "GCF_000888888.8",
    }

    # Results emitted OUT of query order and interleaved.
    # q1 -> 2 hits in genome 7 (counts must increment to 2)
    # q3 -> 1 hit in genome 7, 1 hit in genome 8 (two genome rows)
    # q2 -> 1 hit in genome 8
    # plus a short/garbage line and an unmatched query_id line.
    results = out / "results.full.tsv"
    results.write_text("\n".join([
        fs_line(q3.query_id, "hitY_GCF_000888888.8"),
        fs_line(q1.query_id, "hitX_GCF_000777777.7"),
        "garbage\ttoo\tshort",                                   # D6 skipped
        fs_line(q2.query_id, "hitY_GCF_000888888.8"),
        fs_line(q1.query_id, "hitX_GCF_000777777.7"),
        fs_line("PFZZZ|GCF_x|GHOST|rank9", "hitX_GCF_000777777.7"),  # D5 unmatched
        fs_line(q3.query_id, "hitX_GCF_000777777.7"),
    ]) + "\n")

    write_hits_and_counts(results, queries, target_map, out, progress_interval=0)

    # --- hits.tsv ---
    header, rows = read_tsv(out / "hits.tsv")

    # D2 exact header.
    expected_header = FOLDSEEK_COLUMNS + [
        "target_genome", "pfam_accession", "query_genome",
        "query_target_name", "query_rank", "query_source_file",
    ]
    check(header == expected_header, "D2 hits.tsv header exact and in order")

    # D6 garbage line skipped -> 5 valid hit rows (q3,q1,q2,q1,ghost,q3) minus
    # garbage = 6 written (the ghost line IS valid, just unmatched).
    check(len(rows) == 6, "D6 short line skipped; 6 valid hit rows written")

    by_query = {q.query_id: q for q in queries}
    for row in rows:
        rowmap = dict(zip(header, row))
        # D2 first 12 cells are byte-identical passthrough of the foldseek line.
        # We reconstruct the source line and confirm passthrough.
        passthrough = "\t".join(row[:len(FOLDSEEK_COLUMNS)])
        check(passthrough == fs_line(rowmap["query"], rowmap["target"]),
              f"D2 foldseek columns passed through verbatim for {rowmap['query']}")

        q = by_query.get(rowmap["query"])
        if q is None:
            # D5 unmatched query_id -> NA metadata, no crash.
            check(rowmap["pfam_accession"] == "NA"
                  and rowmap["query_genome"] == "NA"
                  and rowmap["query_target_name"] == "NA"
                  and rowmap["query_rank"] == "NA"
                  and rowmap["query_source_file"] == "NA",
                  "D5 unmatched query_id gets NA metadata")
        else:
            # D1 each row's metadata matches the query for THAT row's query_id.
            check(rowmap["pfam_accession"] == q.pfam_accession
                  and rowmap["query_genome"] == q.query_genome
                  and rowmap["query_target_name"] == q.target_name
                  and rowmap["query_rank"] == q.rank
                  and rowmap["query_source_file"] == q.source_file,
                  f"D1 metadata correctly joined for {q.query_id}")
        # D3 target_genome resolved from the map (distinct from query_genome).
        check(rowmap["target_genome"] == target_map[rowmap["target"]],
              f"D3 target_genome from map for {rowmap['target']}")
        if q is not None:
            check(rowmap["target_genome"] != rowmap["query_genome"],
                  "D3 target_genome and query_genome are not crossed")

    # --- genome_counts.tsv ---
    cheader, crows = read_tsv(out / "genome_counts.tsv")

    # D4 columns == every query_id in original order (incl. zero-hit q4).
    check(cheader == ["genome"] + [q.query_id for q in queries],
          "D4 counts columns == all query_ids in original order (incl zero-hit)")

    counts = {row[0]: dict(zip(cheader[1:], row[1:])) for row in crows}
    # D4 rows sorted by genome.
    check([row[0] for row in crows] == sorted(counts),
          "D4 genome rows sorted by genome name")
    # D4 two genomes present (7 from q1x2+q3, 8 from q2+q3, plus ghost in 7).
    check(set(counts) == {"GCF_000777777.7", "GCF_000888888.8"},
          "D4 exactly the two hit genomes appear as rows")
    # D4 q1 has 2 hits in genome 7.
    check(counts["GCF_000777777.7"][q1.query_id] == "2",
          "D4 repeated (genome,query) hits increment the cell to 2")
    # D4 q3 appears under both genomes (separate rows).
    check(counts["GCF_000777777.7"][q3.query_id] == "1"
          and counts["GCF_000888888.8"][q3.query_id] == "1",
          "D4 same query under two genomes -> two separate rows")
    # D4 zero-hit q4 column all zeros.
    check(all(counts[g][q4.query_id] == "0" for g in counts),
          "D4 zero-hit query column is all zeros")


# --------------------------------------------------------------------------- #
# Group E: end-to-end (gated on foldseek binary)
# --------------------------------------------------------------------------- #
def find_foldseek() -> Path | None:
    env = os.environ.get("FOLDSEEK")
    candidates = [env] if env else []
    candidates.append(str(Path.home() / "foldseek" / "bin" / "foldseek"))
    for c in candidates:
        if c and Path(c).expanduser().exists():
            return Path(c).expanduser()
    return None


def run_pipeline(foldseek: Path, input_tsv: Path, target_db: Path,
                 out_dir: Path, test_flag: bool) -> None:
    cmd = [
        "uv", "run", "python", str(REPO / "foldseek_topn_pfam.py"),
        str(input_tsv), str(target_db), str(out_dir),
        "--foldseek", str(foldseek),
        "--gpu", "0", "--threads", "2", "--max-seqs", "5",
    ]
    if test_flag:
        cmd.append("--test")
    subprocess.run(cmd, check=True, cwd=REPO)


def assert_pipeline_invariants(out_dir: Path) -> None:
    for name in ("hits.tsv", "genome_counts.tsv", "query_metadata.tsv",
                 "results.full.tsv"):
        path = out_dir / name
        check(path.exists() and path.stat().st_size > 0,
              f"E {name} exists and is non-empty")

    mheader, mrows = read_tsv(out_dir / "query_metadata.tsv")
    meta = {row[mheader.index("query_id")]: row for row in mrows}
    pfam_col = mheader.index("pfam_accession")
    genome_col = mheader.index("query_genome")

    hheader, hrows = read_tsv(out_dir / "hits.tsv")
    hq = hheader.index("query")
    hpfam = hheader.index("pfam_accession")
    hgenome = hheader.index("query_genome")
    htgenome = hheader.index("target_genome")

    # Invariant 2: every hit query is in metadata with matching pfam/genome.
    for row in hrows:
        qid = row[hq]
        check(qid in meta, f"E hit query {qid} present in query_metadata.tsv")
        check(row[hpfam] == meta[qid][pfam_col],
              f"E hit pfam matches metadata for {qid}")
        check(row[hgenome] == meta[qid][genome_col],
              f"E hit query_genome matches metadata for {qid}")

    # Invariant 3: genome_counts columns == metadata query_id set.
    cheader, crows = read_tsv(out_dir / "genome_counts.tsv")
    check(set(cheader[1:]) == set(meta),
          "E genome_counts columns == query_id set in metadata")

    # Invariant 4: each genome row's count sum == # hits with that target_genome.
    hits_per_genome: dict[str, int] = {}
    for row in hrows:
        hits_per_genome[row[htgenome]] = hits_per_genome.get(row[htgenome], 0) + 1
    for row in crows:
        genome = row[0]
        row_sum = sum(int(v) for v in row[1:])
        check(row_sum == hits_per_genome.get(genome, 0),
              f"E counts row sum matches hit count for genome {genome}")


def build_real_3di_subset(src: Path, dst: Path, n: int) -> int:
    """Copy the header + first `n` non-empty-3Di rows of `src` into `dst`.

    Gives an end-to-end input guaranteed to contain `n` real full 3Di strings
    while keeping the foldseek search small.
    """
    with src.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        col = header.index("3di_sequence")
        kept = []
        for row in reader:
            s = row[col].strip() if col < len(row) else ""
            if s and s != "NA":
                kept.append(row)
                if len(kept) >= n:
                    break
    with dst.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(kept)
    return len(kept)


def test_e_end_to_end() -> None:
    section("E: end-to-end (real foldseek)")

    foldseek = find_foldseek()
    target_db = REPO / "results_fs_only_full_fs_compare" / "foldseek_db" / "sequenceDB"
    if foldseek is None:
        print("  SKIP: foldseek binary not found")
        return
    if not target_db.exists():
        print(f"  SKIP: target DB not found: {target_db}")
        return

    with tempfile.TemporaryDirectory() as td:
        base = Path(td)

        # A >=10 real full 3Di input. --test would only keep 7 here (the first 3
        # pfams have just 1+3+3 rows with 3Di), so we build a 12-row subset of
        # all_pfam_top10_hits.tsv (query_pfam_accession column variant) instead.
        subset = base / "all_pfam_real12.tsv"
        n_subset = build_real_3di_subset(REPO / "all_pfam_top10_hits.tsv",
                                         subset, 12)
        check(n_subset >= 10, f"E built subset with >=10 real 3Di rows ({n_subset})")

        cases = [
            REPO / "pf00198_top5.tsv",      # domain_accession variant
            REPO / "pfAAA_top10.tsv",       # empty-3Di skipping
            subset,                         # >=10 real 3Di, query_pfam variant
        ]
        for input_tsv in cases:
            if not input_tsv.exists():
                fail(f"E required input file missing: {input_tsv}")
            out_dir = base / (input_tsv.stem + "_out")
            print(f"  running pipeline on {input_tsv.name}")
            run_pipeline(foldseek, input_tsv, target_db, out_dir, test_flag=False)
            assert_pipeline_invariants(out_dir)

            mheader, mrows = read_tsv(out_dir / "query_metadata.tsv")
            pfam_col = mheader.index("pfam_accession")
            seq_col = mheader.index("three_di_sequence")

            if input_tsv is subset:
                # >=10 real full 3Di strings end-to-end.
                check(len(mrows) == n_subset,
                      f"E subset produced all {n_subset} real 3Di queries")
                check(len(mrows) >= 10,
                      f"E subset run covers >=10 real 3Di strings ({len(mrows)})")
                check(all(r[seq_col] and r[seq_col] != "NA" for r in mrows),
                      "E subset: every query has a non-empty 3Di string")
            if input_tsv.name == "pfAAA_top10.tsv":
                # Only 1 of 10 rows has a 3Di sequence; the rest must be dropped.
                check(len(mrows) == 1,
                      "E pfAAA_top10.tsv yields exactly 1 query (9 empty 3Di dropped)")
                check(mrows[0][pfam_col] == "PF00004.36",
                      "E pfAAA_top10.tsv: the single query is PF00004.36")


# --------------------------------------------------------------------------- #
def main() -> int:
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        test_a_load_queries(tmp)
        test_b_query_db(tmp)
        test_c_genome_resolution(tmp)
        test_d_hits_and_counts(tmp)
    test_e_end_to_end()
    print("\nALL TESTS PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
