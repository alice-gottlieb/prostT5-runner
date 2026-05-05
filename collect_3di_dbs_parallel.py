"""
Parallel fork-join version of collect_3di_dbs.py.

Lays merges out as a binary reduction tree: at level 0 we have N task DBs;
at each level we pair them up and run merges concurrently in a thread pool,
feeding results into the next level until one root DB remains.

Concurrency = total_threads // merge_db_threads worker threads, each running
one foldseek concatdbs invocation (which uses merge_db_threads internally).

Fault tolerance: if a pair merge fails, one of its two inputs is propagated
upward as if the merge had not happened (preserves more sequences than
dropping the subtree). The failure is appended to a TSV log. The failed
pair's inputs and any partial output stay on disk for inspection.

Example
-------
    python collect_3di_dbs_parallel.py -b /scratch/all_3dis -o merged_3di_db \\
        -T 16 --merge-db-threads 1
"""

import argparse
import datetime
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from collect_3di_dbs import (
    COMPLETION_MARKER,
    DB_PREFIX,
    DB_SUFFIXES,
    discover_completed_dirs,
    has_full_db,
    merge_metadata,
)


FAILURE_LOG_HEADER = "timestamp\tlevel\tpair_index\tdb_a\tdb_b\tpropagated\terror\n"


# ---------------------------------------------------------------------------
# Pair-merge primitive
# ---------------------------------------------------------------------------

def _db_prefix_in(d: Path) -> str:
    """Original task dirs hold a 'sequenceDB' family; intermediate merge
    outputs hold a 'mergedDB' family. Pick whichever is present."""
    if (d / "mergedDB").exists():
        return "mergedDB"
    if (d / DB_PREFIX).exists():
        return DB_PREFIX
    raise FileNotFoundError(
        f"No '{DB_PREFIX}' or 'mergedDB' family found in {d}")


def merge_pair(a_dir: Path, b_dir: Path, out_dir: Path,
               foldseek_bin: str, threads: int) -> Path:
    """Merge two foldseek DBs into out_dir/mergedDB*. Each input directory
    may hold either a 'sequenceDB' family (original task dirs) or a
    'mergedDB' family (outputs of a previous pair-merge). Raises
    subprocess.CalledProcessError on failure.

    No --preserve-keys (keys must be offset consistently across all three
    sub-DBs to keep main/_h/_ss aligned). No lndb here — the caller does
    one final lndb on the root."""
    out_dir.mkdir(parents=True, exist_ok=True)
    merged = out_dir / "mergedDB"

    a_prefix = _db_prefix_in(a_dir)
    b_prefix = _db_prefix_in(b_dir)

    for src in a_dir.glob(f"{a_prefix}*"):
        dest = out_dir / src.name.replace(a_prefix, "mergedDB", 1)
        shutil.copy2(src, dest)

    next_db = b_dir / b_prefix
    tmp = out_dir / "_concat_tmp"

    cmds = [
        [foldseek_bin, "concatdbs", str(merged), str(next_db),
         str(tmp), "--threads", str(threads)],
        [foldseek_bin, "concatdbs", f"{merged}_h", f"{next_db}_h",
         f"{tmp}_h", "--threads", str(threads)],
        [foldseek_bin, "concatdbs", f"{merged}_ss", f"{next_db}_ss",
         f"{tmp}_ss", "--threads", str(threads)],
    ]
    for cmd in cmds:
        subprocess.run(cmd, check=True, capture_output=True, text=True)

    for old in out_dir.glob("mergedDB*"):
        old.unlink()
    for new in out_dir.glob("_concat_tmp*"):
        new.rename(out_dir / new.name.replace("_concat_tmp", "mergedDB", 1))

    return out_dir


# ---------------------------------------------------------------------------
# Failure log
# ---------------------------------------------------------------------------

def append_failure(log_path: Path, level: int, pair_index: int,
                   a: Path, b: Path, propagated: Path, error: str) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    new_file = not log_path.exists()
    with open(log_path, "a") as f:
        if new_file:
            f.write(FAILURE_LOG_HEADER)
        ts = datetime.datetime.now().isoformat(timespec="seconds")
        err_one_line = error.replace("\t", " ").replace("\n", " | ")[:1000]
        f.write(f"{ts}\t{level}\t{pair_index}\t{a}\t{b}\t{propagated}\t{err_one_line}\n")


# ---------------------------------------------------------------------------
# Tree reduction
# ---------------------------------------------------------------------------

def reduce_tree(leaves: list[Path], out_dir: Path, foldseek_bin: str,
                num_workers: int, merge_threads: int,
                failure_log: Path) -> Path:
    """Iteratively pair-merge `leaves` level-by-level until one DB remains.
    Returns the directory containing the root mergedDB* family."""
    pool = ThreadPoolExecutor(max_workers=num_workers)
    tmp_root = out_dir / "tmp_merge"
    do_not_delete: set[Path] = set()  # any input to a failed merge
    leaf_set = {p.resolve() for p in leaves}

    def safe_cleanup(p: Path) -> None:
        """Delete a temp DB dir, but never the original leaves and never
        anything marked do-not-delete."""
        if p.resolve() in leaf_set:
            return
        if p.resolve() in do_not_delete:
            return
        if not p.is_dir():
            return
        try:
            shutil.rmtree(p)
        except OSError as e:
            print(f"  WARN: could not remove {p}: {e}")

    current = list(leaves)
    level = 1
    while len(current) > 1:
        pairs = [(current[i], current[i + 1])
                 for i in range(0, len(current) - 1, 2)]
        carry = current[-1] if len(current) % 2 == 1 else None

        print(f"\n=== Level {level}: {len(pairs)} pair merges"
              f"{' + 1 carry' if carry else ''} "
              f"({num_workers} workers, {merge_threads} threads each) ===")

        level_root = tmp_root / f"level_{level}"
        level_root.mkdir(parents=True, exist_ok=True)

        futures = []
        for j, (a, b) in enumerate(pairs):
            merge_out = level_root / f"merge_{j}"
            futures.append(pool.submit(
                merge_pair, a, b, merge_out, foldseek_bin, merge_threads))

        next_level: list[Path] = []
        for j, ((a, b), fut) in enumerate(zip(pairs, futures)):
            try:
                result_dir = fut.result()
                next_level.append(result_dir)
                # Inputs successfully consumed — drop them if they were temps.
                safe_cleanup(a)
                safe_cleanup(b)
                print(f"  [L{level} #{j}] OK  {a.name} + {b.name} -> {result_dir}")
            except subprocess.CalledProcessError as e:
                err = (e.stderr or "").strip() or repr(e)
                propagated = a
                do_not_delete.add(a.resolve())
                do_not_delete.add(b.resolve())
                append_failure(failure_log, level, j, a, b, propagated, err)
                next_level.append(propagated)
                print(f"  [L{level} #{j}] FAIL {a.name} + {b.name}: "
                      f"propagating {a.name}; logged to {failure_log}",
                      file=sys.stderr)

        if carry is not None:
            next_level.append(carry)
        current = next_level
        level += 1

    pool.shutdown(wait=True)
    return current[0]


# ---------------------------------------------------------------------------
# Single-DB shortcut and root finalization
# ---------------------------------------------------------------------------

def copy_db_family(src_dir: Path, src_prefix: str,
                   dst_dir: Path, dst_prefix: str) -> None:
    dst_dir.mkdir(parents=True, exist_ok=True)
    for f in src_dir.glob(f"{src_prefix}*"):
        f_name = f.name.replace(src_prefix, dst_prefix, 1)
        shutil.copy2(f, dst_dir / f_name)


def move_db_family(src_dir: Path, src_prefix: str,
                   dst_dir: Path, dst_prefix: str) -> None:
    dst_dir.mkdir(parents=True, exist_ok=True)
    for f in src_dir.glob(f"{src_prefix}*"):
        f_name = f.name.replace(src_prefix, dst_prefix, 1)
        shutil.move(str(f), dst_dir / f_name)


def finalize_root(root_dir: Path, out_dir: Path, foldseek_bin: str) -> Path:
    """Move root tmp_merge mergedDB* files up to out_dir/mergedDB* and run lndb."""
    if root_dir.resolve() != out_dir.resolve():
        # Root lives under tmp_merge — move it up.
        move_db_family(root_dir, "mergedDB", out_dir, "mergedDB")
        try:
            shutil.rmtree(root_dir)
        except OSError:
            pass

    merged = out_dir / "mergedDB"
    subprocess.run(
        [foldseek_bin, "lndb", f"{merged}_h", f"{merged}_ss_h"],
        check=True,
    )

    index_file = out_dir / "mergedDB.index"
    if index_file.exists():
        n = sum(1 for _ in open(index_file))
        print(f"\nMerged DB contains {n} sequences -> {merged}")

    return merged


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("-b", "--base-dir", required=True,
                   help="Directory to scan for completed 3Di outputs")
    p.add_argument("-o", "--output-dir", required=True,
                   help="Directory to write the merged foldseek DB into")
    p.add_argument("-T", "--total-threads", type=int,
                   default=os.cpu_count() or 4,
                   help="Total CPU budget across all concurrent merges "
                        "(default: os.cpu_count())")
    p.add_argument("--merge-db-threads", type=int, default=1,
                   help="Threads per individual pair-merge (default: 1)")
    p.add_argument("--failed-merges-log", default=None,
                   help="Path for the TSV failure log "
                        "(default: <output_dir>/failed_merges.txt)")
    p.add_argument("-f", "--foldseek-path", default="foldseek",
                   help="foldseek binary (default: foldseek on PATH)")
    args = p.parse_args()

    foldseek_bin = args.foldseek_path
    if shutil.which(foldseek_bin) is None and not Path(foldseek_bin).is_file():
        raise SystemExit(f"foldseek not found at '{foldseek_bin}' "
                         f"(use --foldseek-path to override)")

    base = Path(args.base_dir).resolve()
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    failure_log = (Path(args.failed_merges_log) if args.failed_merges_log
                   else out_dir / "failed_merges.txt")

    completed = discover_completed_dirs(
        base)
    if not completed:
        raise SystemExit(
            f"No completed 3Di outputs (no '{COMPLETION_MARKER}') under {base}")

    mergeable, fasta_only = [], []
    for d in completed:
        (mergeable if has_full_db(d) else fasta_only).append(d)

    if fasta_only:
        print(f"\nWARNING: {len(fasta_only)} directories have only "
              f"{COMPLETION_MARKER} but no sibling {DB_PREFIX}* files. Skipping.")
        for d in fasta_only:
            print(f"  - {d}")

    if not mergeable:
        raise SystemExit(
            "\nNo directories with complete sequenceDB files found; nothing to merge.")

    num_workers = max(1, args.total_threads // args.merge_db_threads)
    print(f"\nMerging {len(mergeable)} sequenceDBs into {out_dir} "
          f"({num_workers} workers x {args.merge_db_threads} threads each, "
          f"total budget {args.total_threads})")
    print(f"Failure log: {failure_log}")

    if len(mergeable) == 1:
        print("\nOnly 1 DB found — copying instead of merging.")
        copy_db_family(mergeable[0], DB_PREFIX, out_dir, "mergedDB")
        finalize_root(out_dir, out_dir, foldseek_bin)
    else:
        root = reduce_tree(mergeable, out_dir, foldseek_bin,
                           num_workers, args.merge_db_threads, failure_log)
        finalize_root(root, out_dir, foldseek_bin)
        # Best-effort cleanup of the empty tmp_merge tree.
        tmp_root = out_dir / "tmp_merge"
        if tmp_root.exists():
            try:
                shutil.rmtree(tmp_root)
            except OSError as e:
                print(f"  WARN: could not remove {tmp_root}: {e} "
                      f"(probably contains leftover failed-merge artifacts)")

    merge_metadata(mergeable, out_dir)


if __name__ == "__main__":
    main()
