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
import logging
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from collect_3di_dbs import (
    DB_PREFIX,
    DB_SUFFIXES,
    merge_metadata,
)


# ---------------------------------------------------------------------------
# Discovery (prefix-agnostic — picks up sequenceDB* AND mergedDB* families)
# ---------------------------------------------------------------------------

KNOWN_DB_PREFIXES = (DB_PREFIX, "mergedDB")


def _missing_db_files(d: Path, prefix: str) -> list[str]:
    """Return the list of expected DB-family filenames that don't exist in d."""
    return [f"{prefix}{s}" for s in DB_SUFFIXES
            if not (d / f"{prefix}{s}").exists()]


def has_full_db_any_prefix(d: Path) -> str | None:
    """Return the prefix ('sequenceDB' or 'mergedDB') of a complete foldseek
    DB family in `d`, or None if no complete family is present."""
    for prefix in KNOWN_DB_PREFIXES:
        if not _missing_db_files(d, prefix):
            return prefix
    return None


def discover_db_dirs(base: Path, output_file: Path | None = None) -> list[Path]:
    """Recursively find every directory under `base` that holds a complete
    foldseek DB family of either prefix. Catches per-task outputs
    (`sequenceDB*`) and already-merged outputs (`mergedDB*`).

    Discovery uses os.walk (followlinks=True) so symlinked subdirectories
    are also traversed. Any directory that has a `<prefix>.index` file but
    is missing some other expected suffix is reported with the missing
    filenames so the user can see why it was rejected."""
    base = Path(base).resolve()
    if not base.is_dir():
        raise SystemExit(f"discover_db_dirs: not a directory: {base}")

    candidates: set[Path] = set()
    for root, _, files in os.walk(base, followlinks=True):
        for prefix in KNOWN_DB_PREFIXES:
            if f"{prefix}.index" in files:
                candidates.add(Path(root))
                break  # one match per dir is enough

    accepted: list[Path] = []
    rejected: list[tuple[Path, dict[str, list[str]]]] = []
    for d in sorted(candidates):
        prefix = has_full_db_any_prefix(d)
        if prefix is not None:
            accepted.append(d)
        else:
            missing_per_prefix = {
                p: _missing_db_files(d, p)
                for p in KNOWN_DB_PREFIXES
                if (d / f"{p}.index").exists()
            }
            rejected.append((d, missing_per_prefix))

    print(f"Scanned {base}: {len(candidates)} candidates "
          f"({len(accepted)} complete, {len(rejected)} incomplete) "
          f"[prefixes: {'/'.join(KNOWN_DB_PREFIXES)}]")
    for d in accepted:
        print(f"  OK         {d}")
    for d, missing in rejected:
        for p, miss in missing.items():
            print(f"  INCOMPLETE {d} [{p}] missing: {', '.join(miss)}")

    if output_file is not None:
        output_file.parent.mkdir(parents=True, exist_ok=True)
        with open(output_file, "w") as f:
            for d in accepted:
                f.write(f"{d}\n")
        print(f"Wrote {len(accepted)} discovered DB paths to {output_file}")

    return accepted


FAILURE_LOG_HEADER = "timestamp\tlevel\tpair_index\tdb_a\tdb_b\tpropagated\terror\n"

logger = logging.getLogger("collect_3di_parallel")


def setup_logging(log_file: Path | None) -> None:
    """Configure the module logger with a stream handler (stderr) and an
    optional file handler. The logging module's handlers are thread-safe
    (each emit() acquires an internal RLock), so worker threads can call
    logger.info(...) concurrently without garbled output."""
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()
    fmt = logging.Formatter(
        "%(asctime)s [%(threadName)s] %(levelname)-7s %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
    )

    stream_h = logging.StreamHandler(sys.stderr)
    stream_h.setFormatter(fmt)
    stream_h.setLevel(logging.INFO)
    logger.addHandler(stream_h)

    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        file_h = logging.FileHandler(log_file, mode="a")
        file_h.setFormatter(fmt)
        file_h.setLevel(logging.DEBUG)
        logger.addHandler(file_h)
        logger.info("Logging to %s", log_file)


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

def _threads_per_merge(total_threads: int, n_merges: int,
                       min_per_merge: int) -> list[int]:
    """Distribute total_threads across n_merges. Each merge gets at least
    min_per_merge. If total_threads >= n_merges * min_per_merge, the spare
    budget is spread across merges (some get one extra). Otherwise we cannot
    fit all merges with the floor, so each gets exactly min_per_merge and
    the pool will throttle concurrency."""
    if n_merges == 0:
        return []
    base = max(min_per_merge, total_threads // n_merges)
    if base * n_merges <= total_threads:
        threads = [base] * n_merges
        remainder = total_threads - base * n_merges
        for i in range(remainder):
            threads[i] += 1
        return threads
    return [min_per_merge] * n_merges


def reduce_tree(leaves: list[Path], out_dir: Path, foldseek_bin: str,
                total_threads: int, min_merge_threads: int,
                failure_log: Path) -> Path:
    """Iteratively pair-merge `leaves` level-by-level until one DB remains.
    Returns the directory containing the root mergedDB* family.

    Threads per merge are computed dynamically per level: if the level's
    pair count <= total_threads, the budget is split across merges (so a
    near-final level with few merges gets fat-thread allocations); if not,
    each merge takes the floor and the pool throttles concurrency."""
    pool = ThreadPoolExecutor(max_workers=max(1, total_threads))
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
            logger.warning("could not remove %s: %s", p, e)

    current = list(leaves)
    level = 1
    while len(current) > 1:
        pairs = [(current[i], current[i + 1])
                 for i in range(0, len(current) - 1, 2)]
        carry = current[-1] if len(current) % 2 == 1 else None

        thread_alloc = _threads_per_merge(total_threads, len(pairs),
                                          min_merge_threads)
        if thread_alloc and len(set(thread_alloc)) == 1:
            alloc_desc = f"{thread_alloc[0]} threads each"
        elif thread_alloc:
            alloc_desc = (f"{min(thread_alloc)}-{max(thread_alloc)} threads each "
                          f"(sum={sum(thread_alloc)})")
        else:
            alloc_desc = "no merges"
        logger.info("=== Level %d: %d pair merges%s (budget=%d, %s) ===",
                    level, len(pairs),
                    " + 1 carry" if carry else "",
                    total_threads, alloc_desc)

        level_root = tmp_root / f"level_{level}"
        level_root.mkdir(parents=True, exist_ok=True)

        futures = []
        for j, (a, b) in enumerate(pairs):
            merge_out = level_root / f"merge_{j}"
            futures.append(pool.submit(
                merge_pair, a, b, merge_out, foldseek_bin, thread_alloc[j]))

        next_level: list[Path] = []
        for j, ((a, b), fut) in enumerate(zip(pairs, futures)):
            try:
                result_dir = fut.result()
                next_level.append(result_dir)
                # Inputs successfully consumed — drop them if they were temps.
                safe_cleanup(a)
                safe_cleanup(b)
                logger.info("[L%d #%d] OK  %s + %s -> %s",
                            level, j, a.name, b.name, result_dir)
            except subprocess.CalledProcessError as e:
                err = (e.stderr or "").strip() or repr(e)
                propagated = a
                do_not_delete.add(a.resolve())
                do_not_delete.add(b.resolve())
                append_failure(failure_log, level, j, a, b, propagated, err)
                next_level.append(propagated)
                logger.error("[L%d #%d] FAIL %s + %s: propagating %s; logged to %s",
                             level, j, a.name, b.name, a.name, failure_log)

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
        logger.info("Merged DB contains %d sequences -> %s", n, merged)

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
                   help="Minimum threads per individual pair-merge "
                        "(default: 1). When fewer merges than --total-threads "
                        "remain at a level, the spare budget is divided "
                        "across them so each gets more than this floor.")
    p.add_argument("--failed-merges-log", default=None,
                   help="Path for the TSV failure log "
                        "(default: <output_dir>/failed_merges.txt)")
    p.add_argument("--log-file", default=None,
                   help="Path for the run log "
                        "(default: <output_dir>/run.log). "
                        "Pass an empty string to disable file logging.")
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

    if args.log_file is None:
        log_file = out_dir / "run.log"
    elif args.log_file == "":
        log_file = None
    else:
        log_file = Path(args.log_file)
    setup_logging(log_file)

    mergeable = discover_db_dirs(base, output_file=out_dir / "discovered_dbs.txt")
    if not mergeable:
        raise SystemExit(
            f"No directories with a complete foldseek DB family "
            f"({'/'.join(KNOWN_DB_PREFIXES)}) under {base}")

    logger.info("Merging %d sequenceDBs into %s (total budget %d, min %d threads/merge)",
                len(mergeable), out_dir, args.total_threads, args.merge_db_threads)
    logger.info("Failure log: %s", failure_log)

    if len(mergeable) == 1:
        logger.info("Only 1 DB found — copying instead of merging.")
        src_prefix = has_full_db_any_prefix(mergeable[0]) or DB_PREFIX
        copy_db_family(mergeable[0], src_prefix, out_dir, "mergedDB")
        finalize_root(out_dir, out_dir, foldseek_bin)
    else:
        root = reduce_tree(mergeable, out_dir, foldseek_bin,
                           args.total_threads, args.merge_db_threads,
                           failure_log)
        finalize_root(root, out_dir, foldseek_bin)
        # Best-effort cleanup of the empty tmp_merge tree.
        tmp_root = out_dir / "tmp_merge"
        if tmp_root.exists():
            try:
                shutil.rmtree(tmp_root)
            except OSError as e:
                logger.warning("could not remove %s: %s "
                               "(probably contains leftover failed-merge artifacts)",
                               tmp_root, e)

    merge_metadata(mergeable, out_dir)


if __name__ == "__main__":
    main()
