"""
Smoke test for collect_3di_dbs.py / merge_two_runs.py.

Runs the full merge pipeline against two finished submit_3di_array.sh runs
and verifies:
  1. all DBs are discovered correctly,
  2. originals are not modified during the merge,
  3. the merged DB has the expected sequence count,
  4. the merged 3Di FASTA contains the same headers and sequences as the
     per-task all_sequences_3di.fasta files,
  5. the merged metadata.json was written and is consistent.

Exits non-zero on the first failed check.

Usage:
    python test_merge.py <run_a> <run_b> <output_dir> [--foldseek <bin>]
"""

import argparse
import glob
import hashlib
import json
import subprocess
import sys
from pathlib import Path

from Bio import SeqIO

from collect_3di_dbs import discover_completed_dirs, has_full_db


def section(title: str) -> None:
    print(f"\n=== {title} ===")


def fail(msg: str) -> None:
    print(f"FAIL: {msg}", file=sys.stderr)
    sys.exit(1)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def hash_source_dbs(runs: list[Path]) -> dict[Path, str]:
    """Hash every sequenceDB* file under any foldseek_db/ in the source runs."""
    hashes: dict[Path, str] = {}
    for run in runs:
        for db_dir in run.rglob("foldseek_db"):
            for f in db_dir.glob("sequenceDB*"):
                if f.is_file():
                    hashes[f] = sha256(f)
    return hashes


def load_3di_fasta(paths: list[str]) -> dict[str, str]:
    """Parse {id: 3Di sequence} from a list of FASTA paths."""
    out: dict[str, str] = {}
    for p in paths:
        for rec in SeqIO.parse(p, "fasta"):
            out[rec.id] = str(rec.seq)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_a", type=Path)
    ap.add_argument("run_b", type=Path)
    ap.add_argument("output", type=Path,
                    help="Output dir for the merged DB (must not already contain a merged DB)")
    ap.add_argument("--foldseek", default="foldseek",
                    help="foldseek binary (default: foldseek on PATH)")
    args = ap.parse_args()

    repo_dir = Path(__file__).resolve().parent

    # 1) Discover all DBs
    section("1) Discovering DBs")
    for run in (args.run_a, args.run_b):
        dirs = discover_completed_dirs(run)
        print(f"\n{run}: {len(dirs)} completed")
        for d in dirs:
            tag = "OK         " if has_full_db(d) else "FASTA-only "
            print(f"  {tag}{d}")

    # 2) Snapshot originals, then run the merge
    section("2) Hashing originals + running merge")
    orig_hashes = hash_source_dbs([args.run_a, args.run_b])
    print(f"Hashed {len(orig_hashes)} source files")

    subprocess.run(
        [sys.executable, str(repo_dir / "merge_two_runs.py"),
         str(args.run_a), str(args.run_b), str(args.output),
         "-f", args.foldseek],
        check=True,
    )

    # 3) Confirm originals were untouched
    section("3) Verifying originals unmodified")
    changed = [p for p, h in orig_hashes.items() if not p.exists() or sha256(p) != h]
    if changed:
        for p in changed:
            print(f"  CHANGED: {p}", file=sys.stderr)
        fail("at least one original file changed")
    print("All originals intact ✓")

    # 4a) Index line counts
    section("4a) Index line counts")
    src_count = 0
    for run in (args.run_a, args.run_b):
        for idx in run.rglob("foldseek_db/sequenceDB.index"):
            with open(idx) as f:
                src_count += sum(1 for _ in f)
    with open(args.output / "mergedDB.index") as f:
        merged_count = sum(1 for _ in f)
    print(f"source: {src_count}   merged: {merged_count}")
    if src_count != merged_count:
        fail("merged index count does not match source total")
    print("Counts match ✓")

    # 4b) 3Di FASTA header diff
    section("4b) 3Di FASTA header diff")
    merged_3di = args.output / "merged_3di.fasta"
    subprocess.run(
        [args.foldseek, "convert2fasta",
         str(args.output / "mergedDB_ss"), str(merged_3di)],
        check=True,
    )
    merged_headers = {r.id for r in SeqIO.parse(merged_3di, "fasta")}

    source_fastas: list[str] = []
    for run in (args.run_a, args.run_b):
        source_fastas.extend(str(p) for p in run.rglob("all_sequences_3di.fasta"))
    source_headers = {r.id for f in source_fastas for r in SeqIO.parse(f, "fasta")}

    only_src = sorted(source_headers - merged_headers)
    only_merged = sorted(merged_headers - source_headers)
    if only_src or only_merged:
        for h in only_src[:5]:
            print(f"  in source only: {h}")
        for h in only_merged[:5]:
            print(f"  in merged only: {h}")
        fail(f"header sets differ ({len(only_src)} src-only, {len(only_merged)} merged-only)")
    print(f"Headers match ✓ ({len(merged_headers)} sequences)")

    # 4c) 3Di sequence round-trip
    section("4c) 3Di sequence round-trip")
    src_seqs = load_3di_fasta(source_fastas)
    mismatches = 0
    for r in SeqIO.parse(merged_3di, "fasta"):
        if src_seqs.get(r.id) != str(r.seq):
            mismatches += 1
            if mismatches <= 5:
                print(f"  MISMATCH {r.id}")
    print(f"{mismatches} mismatches across {len(src_seqs)} sequences")
    if mismatches:
        fail("3Di sequences do not round-trip cleanly")

    # 4d) Merged metadata.json
    section("4d) Merged metadata.json")
    meta_path = args.output / "metadata.json"
    if not meta_path.exists():
        fail(f"{meta_path} was not written")
    with open(meta_path) as f:
        meta = json.load(f)
    print(f"source runs: {len(meta['source_runs'])}")
    print(f"unique downloaded assemblies: {len(meta['assemblies_downloaded'])}")
    print(f"total num_sequences: {meta['num_sequences']}")

    section("All checks passed ✓")


if __name__ == "__main__":
    main()
