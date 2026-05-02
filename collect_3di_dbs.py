"""
Find all completed 3Di outputs under a directory and concatenate them into a
single foldseek database using `foldseek concatdbs`.

A "completed" run is identified by the presence of `all_sequences_3di.fasta`
inside a foldseek_db-style directory (the file that batch_3di_foldseek.py
writes after it finishes extracting 3Di codes).

Why we merge sequenceDB files instead of the 3Di FASTAs themselves
------------------------------------------------------------------
`foldseek createdb` accepts only:
    - amino acid FASTA (then runs ProstT5 to predict 3Di), or
    - structure files (PDB / mmCIF, then derives 3Di from coordinates).

It does NOT accept a pre-computed 3Di FASTA as input. A foldseek searchable
DB requires the paired AA DB, 3Di DB (`_ss`), and headers DB (`_h`) — the
3Di FASTA on its own is a lossy export. So if a discovered run has only
`all_sequences_3di.fasta` and no sibling `sequenceDB*` files, this script
cannot fold it into the merged DB and will warn instead.

When sequenceDB files are present (the normal case for batch_3di_foldseek.py
output), `foldseek concatdbs` (an MMseqs2 module foldseek re-exports) is
used to iteratively merge the main, header, and 3Di sub-databases.

Example
-------
    python collect_3di_dbs.py -b /scratch/all_3dis -o merged_3di_db
"""

import argparse
import json
import os
import shutil
import subprocess
from pathlib import Path


COMPLETION_MARKER = "all_sequences_3di.fasta"
DB_PREFIX = "sequenceDB"
DB_SUFFIXES = [
    "", ".dbtype", ".index",
    "_h", "_h.dbtype", "_h.index",
    "_ss", "_ss.dbtype", "_ss.index",
]


def run(cmd: list[str], description: str) -> None:
    print(f"\n{description}: {' '.join(cmd)}", flush=True)
    subprocess.run(cmd, check=True)


def discover_completed_dirs(base: Path) -> list[Path]:
    """Find every directory containing a completed 3Di FASTA."""
    dirs = []
    for root, _, files in os.walk(base):
        if COMPLETION_MARKER in files:
            dirs.append(Path(root))
    dirs.sort()
    return dirs


def has_full_db(db_dir: Path) -> bool:
    """True if all foldseek sequenceDB files for merging are present."""
    return all((db_dir / f"{DB_PREFIX}{s}").exists() for s in DB_SUFFIXES)


def merge_metadata(db_dirs: list[Path], out_dir: Path) -> Path | None:
    """Aggregate metadata.json from each task dir (parent of foldseek_db/)
    into a single merged metadata.json in out_dir. Accession lists are
    deduplicated; numeric fields are summed."""
    merged = {
        "source_runs": [],
        "assemblies_requested": [],
        "assemblies_downloaded": [],
        "num_sequences": 0,
    }
    seen_requested: set[str] = set()
    seen_downloaded: set[str] = set()
    found_any = False

    for db_dir in db_dirs:
        meta_path = db_dir.parent / "metadata.json"
        if not meta_path.exists():
            print(f"  WARN: no metadata.json next to {db_dir}, skipping")
            continue
        try:
            with open(meta_path) as f:
                meta = json.load(f)
        except (json.JSONDecodeError, OSError) as e:
            print(f"  WARN: could not read {meta_path}: {e}")
            continue

        found_any = True
        merged["source_runs"].append(str(meta_path))
        for acc in meta.get("assemblies_requested", []):
            if acc not in seen_requested:
                seen_requested.add(acc)
                merged["assemblies_requested"].append(acc)
        for acc in meta.get("assemblies_downloaded", []):
            if acc not in seen_downloaded:
                seen_downloaded.add(acc)
                merged["assemblies_downloaded"].append(acc)
        merged["num_sequences"] += int(meta.get("num_sequences") or 0)

    if not found_any:
        return None

    out_path = out_dir / "metadata.json"
    with open(out_path, "w") as f:
        json.dump(merged, f, indent=2)
    print(f"\nMerged metadata: {len(merged['assemblies_downloaded'])} unique "
          f"assemblies across {len(merged['source_runs'])} runs -> {out_path}")
    return out_path


def merge_dbs(db_dirs: list[Path], out_dir: Path, foldseek_bin: str, threads: int) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    merged = out_dir / "mergedDB"

    # Seed merged DB with the first run's files.
    first = db_dirs[0]
    print(f"Seeding merged DB from {first}")
    for src in first.glob(f"{DB_PREFIX}*"):
        dest = out_dir / src.name.replace(DB_PREFIX, "mergedDB", 1)
        shutil.copy2(src, dest)

    tmp = out_dir / "_concat_tmp"

    for i, db_dir in enumerate(db_dirs[1:], start=2):
        next_db = db_dir / DB_PREFIX
        print(f"\n[{i}/{len(db_dirs)}] Concatenating {db_dir}")

        run([foldseek_bin, "concatdbs", str(merged), str(next_db),
             str(tmp), "--threads", str(threads)],
            "  main DB")
        run([foldseek_bin, "concatdbs", f"{merged}_h", f"{next_db}_h",
             f"{tmp}_h", "--threads", str(threads)],
            "  header DB")
        run([foldseek_bin, "concatdbs", f"{merged}_ss", f"{next_db}_ss",
             f"{tmp}_ss", "--threads", str(threads)],
            "  3Di DB")

        # Replace merged DB files with the freshly-concatenated ones.
        for old in out_dir.glob("mergedDB*"):
            old.unlink()
        for new in out_dir.glob("_concat_tmp*"):
            new.rename(out_dir / new.name.replace("_concat_tmp", "mergedDB", 1))

    # Re-link headers onto the 3Di sub-DB so convert2fasta works.
    run([foldseek_bin, "lndb", f"{merged}_h", f"{merged}_ss_h"],
        "Linking headers to 3Di sub-DB")

    index_file = out_dir / "mergedDB.index"
    if index_file.exists():
        n = sum(1 for _ in open(index_file))
        print(f"\nMerged DB contains {n} sequences -> {merged}")

    return merged


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-b", "--base-dir", required=True,
                   help="Directory to scan for completed 3Di outputs")
    p.add_argument("-o", "--output-dir", required=True,
                   help="Directory to write the merged foldseek DB into")
    p.add_argument("-f", "--foldseek-path", default="foldseek",
                   help="foldseek binary (default: foldseek on PATH)")
    p.add_argument("-t", "--threads", type=int, default=os.cpu_count() or 4,
                   help="Threads for concatdbs (default: all CPUs)")
    args = p.parse_args()

    foldseek_bin = args.foldseek_path
    if shutil.which(foldseek_bin) is None and not Path(foldseek_bin).is_file():
        raise SystemExit(f"foldseek not found at '{foldseek_bin}' "
                         f"(use --foldseek-path to override)")

    base = Path(args.base_dir).resolve()
    completed = discover_completed_dirs(base)
    if not completed:
        raise SystemExit(f"No completed 3Di outputs (no '{COMPLETION_MARKER}') under {base}")

    print(f"Found {len(completed)} completed 3Di output directories under {base}")

    mergeable, fasta_only = [], []
    for d in completed:
        (mergeable if has_full_db(d) else fasta_only).append(d)

    if fasta_only:
        print(f"\nWARNING: {len(fasta_only)} directories have only "
              f"{COMPLETION_MARKER} but no sibling {DB_PREFIX}* files.")
        print("These cannot be folded into a foldseek DB: foldseek createdb "
              "does not accept 3Di-only FASTA as input (it requires AA FASTA "
              "or structures). Re-run batch_3di_foldseek.py for these to "
              "regenerate the sequenceDB.")
        for d in fasta_only:
            print(f"  - {d}")

    if not mergeable:
        raise SystemExit("\nNo directories with complete sequenceDB files found; nothing to merge.")

    print(f"\nMerging {len(mergeable)} sequenceDBs into {args.output_dir}")
    out_dir = Path(args.output_dir)
    merge_dbs(mergeable, out_dir, foldseek_bin, args.threads)
    merge_metadata(mergeable, out_dir)


if __name__ == "__main__":
    main()
