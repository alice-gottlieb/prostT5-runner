"""
Merge two finished `submit_3di_array.sh` runs into a single foldseek DB.

Each input folder is the OUTPUT_BASE of one array job — i.e. it contains
`task_NN/foldseek_db/sequenceDB*` directories (one per array task). This
script discovers every completed task DB under both folders and concatenates
them into one foldseek DB at the third path.

Tarballs are not extracted — if a task is only on disk as `task_NN.tar.gz`,
extract it first (e.g. `tar -xzf task_NN.tar.gz`).

Example
-------
    python merge_two_runs.py /scratch/run_a /scratch/run_b /scratch/merged
"""

import argparse
import os
from pathlib import Path

from collect_3di_dbs import discover_completed_dirs, has_full_db, merge_dbs


def collect_task_dbs(run_dir: Path) -> tuple[list[Path], list[Path]]:
    """Return (mergeable_db_dirs, fasta_only_dirs) found under a run folder."""
    completed = discover_completed_dirs(run_dir)
    mergeable, fasta_only = [], []
    for d in completed:
        (mergeable if has_full_db(d) else fasta_only).append(d)
    return mergeable, fasta_only


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("run_a", type=Path, help="First finished run folder (OUTPUT_BASE)")
    p.add_argument("run_b", type=Path, help="Second finished run folder (OUTPUT_BASE)")
    p.add_argument("output", type=Path, help="Output folder for the merged DB")
    p.add_argument("-f", "--foldseek-path", default="foldseek",
                   help="foldseek binary (default: foldseek on PATH)")
    p.add_argument("-t", "--threads", type=int, default=os.cpu_count() or 4,
                   help="Threads for concatdbs (default: all CPUs)")
    args = p.parse_args()

    for run in (args.run_a, args.run_b):
        if not run.is_dir():
            raise SystemExit(f"Not a directory: {run}")

    mergeable_a, fasta_only_a = collect_task_dbs(args.run_a)
    mergeable_b, fasta_only_b = collect_task_dbs(args.run_b)

    print(f"Run A ({args.run_a}): {len(mergeable_a)} mergeable, "
          f"{len(fasta_only_a)} FASTA-only (skipped)")
    print(f"Run B ({args.run_b}): {len(mergeable_b)} mergeable, "
          f"{len(fasta_only_b)} FASTA-only (skipped)")

    for d in fasta_only_a + fasta_only_b:
        print(f"  WARN: skipping {d} (no sequenceDB family)")

    all_dbs = mergeable_a + mergeable_b
    if not all_dbs:
        raise SystemExit("No mergeable task DBs found in either run folder.")

    print(f"\nMerging {len(all_dbs)} task DBs into {args.output}")
    merge_dbs(all_dbs, args.output, args.foldseek_path, args.threads)


if __name__ == "__main__":
    main()
