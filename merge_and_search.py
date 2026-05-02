"""
Merge multiple foldseek databases (from task directories containing metadata.json)
into a single database and run an all-vs-all structural search.

Input:  task directories produced by batch_3di_foldseek.py, each containing
        metadata.json and a foldseek_db/ subdirectory with sequenceDB files.
Output: a merged foldseek database and all-vs-all search results in TSV format.

Requires: foldseek (v9+ with ProstT5 support)

Example usages:

    # Auto-find all task dirs under a base directory
    python merge_and_search.py -o merged_output -b /path/to/results

    # Specify directories explicitly
    python merge_and_search.py -o merged_output -d task_001 task_002 task_003

    # Merge only (skip search)
    python merge_and_search.py -o merged_output -b /path/to/results --merge-only

    # Search only (on already-merged DB)
    python merge_and_search.py -o merged_output --search-only

    # Custom search parameters
    python merge_and_search.py -o merged_output -b /path/to/results \\
        -e 1e-5 -s 7.5 -t 32 --max-seqs 2000

    # Pass extra arguments to foldseek search
    python merge_and_search.py -o merged_output -b /path/to/results \\
        --search-args --alignment-type 2
"""

import argparse
import json
import os
import shutil
import subprocess
import time
from pathlib import Path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _run(cmd: list[str], description: str, stream: bool = False):
    """Run a subprocess command with logging."""
    print(f"\n{description}: {' '.join(cmd)}", flush=True)
    try:
        if stream:
            subprocess.run(cmd, check=True)
        else:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            if result.stdout:
                print(result.stdout[:1000])
            if result.stderr:
                print("Info:", result.stderr[:500])
    except subprocess.CalledProcessError as e:
        if not stream:
            print(f"Command failed:\n  stdout: {e.stdout}\n  stderr: {e.stderr}")
        raise


def get_default_threads() -> int:
    """Return number of available CPUs."""
    try:
        return os.cpu_count() or 4
    except Exception:
        return 4


# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------

DB_SUFFIXES = ["", ".dbtype", ".index",
               "_h", "_h.dbtype", "_h.index",
               "_ss", "_ss.dbtype", "_ss.index"]


def discover_task_dirs(base_dir: str, max_depth: int = 3) -> list[Path]:
    """Find all directories under base_dir that contain metadata.json."""
    base = Path(base_dir).resolve()
    dirs = []
    for root, _, files in os.walk(base):
        # Enforce max depth
        depth = len(Path(root).relative_to(base).parts)
        if depth > max_depth:
            continue
        if "metadata.json" in files:
            dirs.append(Path(root))
    dirs.sort()
    return dirs


def validate_db(task_dir: Path, db_subdir: str, db_prefix: str) -> bool:
    """Check that a task directory has all required foldseek DB files."""
    db_path = task_dir / db_subdir / db_prefix
    for suffix in DB_SUFFIXES:
        if not (db_path.parent / f"{db_path.name}{suffix}").exists():
            print(f"  Warning: missing {db_path}{suffix}, skipping {task_dir.name}")
            return False
    return True


# ---------------------------------------------------------------------------
# Merge
# ---------------------------------------------------------------------------

def merge_databases(
    task_dirs: list[Path],
    merged_dir: Path,
    foldseek_bin: str,
    threads: int,
    db_subdir: str = "foldseek_db",
    db_prefix: str = "sequenceDB",
) -> Path:
    """
    Iteratively merge foldseek databases from multiple task directories.

    Returns the path prefix of the merged database.
    """
    # Validate
    valid_dirs = [d for d in task_dirs if validate_db(d, db_subdir, db_prefix)]
    if not valid_dirs:
        raise SystemExit("Error: no valid databases found")

    merged_dir.mkdir(parents=True, exist_ok=True)
    merged = merged_dir / "mergedDB"

    if len(valid_dirs) == 1:
        print(f"Only 1 valid DB found ({valid_dirs[0].name}). Copying instead of merging.")
        src = valid_dirs[0] / db_subdir / db_prefix
        for f in src.parent.glob(f"{db_prefix}*"):
            dest = merged_dir / f.name.replace(db_prefix, "mergedDB", 1)
            shutil.copy2(f, dest)
        return merged

    print(f"Merging {len(valid_dirs)} databases...")

    # Copy first DB as starting point
    first = valid_dirs[0] / db_subdir / db_prefix
    for f in first.parent.glob(f"{db_prefix}*"):
        dest = merged_dir / f.name.replace(db_prefix, "mergedDB", 1)
        shutil.copy2(f, dest)

    tmp_a = merged_dir / "_tmp_a"

    for i, task_dir in enumerate(valid_dirs[1:], start=2):
        next_db = task_dir / db_subdir / db_prefix
        print(f"[{i}/{len(valid_dirs)}] Merging {task_dir.name}...")

        # Merge main sequence DB
        _run([foldseek_bin, "concatdbs",
              str(merged), str(next_db), str(tmp_a),
              "--threads", str(threads)],
             f"  Merging main DB")

        # Merge header DB
        _run([foldseek_bin, "concatdbs",
              f"{merged}_h", f"{next_db}_h", f"{tmp_a}_h",
              "--threads", str(threads)],
             f"  Merging header DB")

        # Merge 3Di DB
        _run([foldseek_bin, "concatdbs",
              f"{merged}_ss", f"{next_db}_ss", f"{tmp_a}_ss",
              "--threads", str(threads)],
             f"  Merging 3Di DB")

        # Swap: remove old merged, rename tmp to merged
        for suffix in DB_SUFFIXES + ["_ss_h", "_ss_h.dbtype", "_ss_h.index",
                                      ".lookup", ".source"]:
            path = merged.parent / f"mergedDB{suffix}"
            if path.exists():
                path.unlink()

        for f in merged_dir.glob("_tmp_a*"):
            dest = merged_dir / f.name.replace("_tmp_a", "mergedDB", 1)
            f.rename(dest)

    # Re-link 3Di headers
    _run([foldseek_bin, "lndb", f"{merged}_h", f"{merged}_ss_h"],
         "Linking 3Di headers")

    # Count sequences
    index_file = merged.parent / "mergedDB.index"
    if index_file.exists():
        num_seqs = sum(1 for _ in open(index_file))
        print(f"Merge complete: {num_seqs} sequences in merged database")

    return merged


# ---------------------------------------------------------------------------
# Search
# ---------------------------------------------------------------------------

def run_search(
    merged: Path,
    results_dir: Path,
    foldseek_bin: str,
    threads: int,
    evalue: str = "0.001",
    max_seqs: int = 1000,
    sensitivity: float | None = None,
    search_args: list[str] | None = None,
) -> Path:
    """
    Run foldseek all-vs-all search on a merged database and convert results
    to TSV format.

    Returns the path to the results TSV.
    """
    if not merged.exists():
        raise SystemExit(f"Error: merged DB not found at {merged}")

    results_dir.mkdir(parents=True, exist_ok=True)
    result_prefix = results_dir / "result"
    result_tsv = results_dir / "all_vs_all_results.tsv"
    tmp_dir = results_dir / "tmp_search"
    tmp_dir.mkdir(exist_ok=True)

    # Build search command
    cmd = [
        foldseek_bin, "search",
        str(merged), str(merged),
        str(result_prefix), str(tmp_dir),
        "--threads", str(threads),
        "-e", evalue,
        "--max-seqs", str(max_seqs),
    ]
    if sensitivity is not None:
        cmd.extend(["-s", str(sensitivity)])
    if search_args:
        cmd.extend(search_args)

    start = time.time()
    _run(cmd, f"Running all-vs-all search with {threads} threads", stream=True)
    elapsed = time.time() - start
    print(f"Search completed in {elapsed:.0f}s")

    # Convert to TSV (standard format)
    _run([foldseek_bin, "convertalis",
          str(merged), str(merged),
          str(result_prefix), str(result_tsv),
          "--threads", str(threads)],
         "Converting results to TSV")

    num_hits = sum(1 for _ in open(result_tsv))
    print(f"Results: {num_hits} hits written to {result_tsv}")

    return result_tsv


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Merge foldseek databases and run all-vs-all search",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Output
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output directory for merged DB and results",
    )

    # Input (mutually exclusive)
    input_group = parser.add_mutually_exclusive_group()
    input_group.add_argument(
        "-d", "--dirs", nargs="+",
        help="Task directories to merge (each must contain metadata.json)",
    )
    input_group.add_argument(
        "-b", "--base-dir", default=".",
        help="Base directory to scan for subdirectories containing metadata.json "
             "(default: current directory)",
    )

    # Foldseek options
    parser.add_argument(
        "-f", "--foldseek-path", default="foldseek",
        help="Path to foldseek binary (default: foldseek)",
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=get_default_threads(),
        help=f"Threads for search (default: {get_default_threads()})",
    )
    parser.add_argument(
        "-e", "--evalue", default="0.001",
        help="E-value threshold for search (default: 0.001)",
    )
    parser.add_argument(
        "--max-seqs", type=int, default=1000,
        help="Max hits per query (default: 1000)",
    )
    parser.add_argument(
        "-s", "--sensitivity", type=float, default=None,
        help="Sensitivity 1.0-9.5 (default: foldseek default ~5.7)",
    )
    parser.add_argument(
        "--search-args", nargs=argparse.REMAINDER, default=[],
        help="Extra arguments forwarded to foldseek search",
    )

    # Workflow control
    parser.add_argument(
        "--merge-only", action="store_true",
        help="Merge databases only, skip search",
    )
    parser.add_argument(
        "--search-only", action="store_true",
        help="Skip merge, run search on existing merged DB in output dir",
    )

    # DB layout
    parser.add_argument(
        "--db-prefix", default="sequenceDB",
        help="DB prefix inside each task dir (default: sequenceDB)",
    )
    parser.add_argument(
        "--db-subdir", default="foldseek_db",
        help="Subdirectory containing the DB inside each task dir (default: foldseek_db)",
    )

    args = parser.parse_args()

    out = Path(args.output)
    merged_dir = out / "merged_db"
    results_dir = out / "search_results"
    foldseek_bin = args.foldseek_path

    # Check foldseek is available
    if shutil.which(foldseek_bin) is None and not os.path.isfile(foldseek_bin):
        raise SystemExit(
            f"Error: foldseek not found at '{foldseek_bin}'\n"
            "Specify path with --foldseek-path /path/to/foldseek"
        )

    out.mkdir(parents=True, exist_ok=True)

    print("=" * 50)
    print("  Foldseek DB Merge & All-vs-All Search")
    print("=" * 50)

    if not args.search_only:
        # Discover task directories
        if args.dirs:
            task_dirs = []
            for d in args.dirs:
                p = Path(d).resolve()
                if (p / "metadata.json").exists():
                    task_dirs.append(p)
                else:
                    print(f"Warning: {d} has no metadata.json, skipping")
        else:
            task_dirs = discover_task_dirs(args.base_dir)

        if not task_dirs:
            raise SystemExit("Error: no directories with metadata.json found")

        print(f"Found {len(task_dirs)} task directories with metadata.json")

        merged = merge_databases(
            task_dirs, merged_dir, foldseek_bin, args.threads,
            db_subdir=args.db_subdir, db_prefix=args.db_prefix,
        )
    else:
        print("Skipping merge (--search-only)")
        merged = merged_dir / "mergedDB"

    if not args.merge_only:
        result_tsv = run_search(
            merged, results_dir, foldseek_bin, args.threads,
            evalue=args.evalue, max_seqs=args.max_seqs,
            sensitivity=args.sensitivity,
            search_args=args.search_args or None,
        )

    print()
    print("=" * 50)
    print("  Done!")
    print(f"  Merged DB:  {merged_dir}/mergedDB")
    if not args.merge_only:
        print(f"  Results:    {results_dir}/all_vs_all_results.tsv")
    print("=" * 50)


if __name__ == "__main__":
    main()
