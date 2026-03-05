"""
Benchmark batch_3di_foldseek.py: measure wall-clock time per genome across
different thread counts.

Produces a CSV with columns:
    accession, num_proteins, total_residues, threads, time_createdb_s, time_search_s, time_total_s

Example usages:

    # Benchmark all assemblies in the file with thread counts 1,2,4,8
    uv run python benchmark_foldseek.py assemblies.txt \
        --foldseek-path foldseek/bin/foldseek \
        --prostt5-weights ~/prostT5-runner/prostt5_weights/prostt5/prostt5-f16.gguf \
        --gpu --thread-counts 1 2 4 8

    # Benchmark with more thread counts and 3 repetitions for statistics
    uv run python benchmark_foldseek.py assemblies.txt \
        --foldseek-path foldseek/bin/foldseek \
        --prostt5-weights ~/prostT5-runner/prostt5_weights/prostt5/prostt5-f16.gguf \
        --gpu --thread-counts 1 2 4 8 16 --reps 3

    # Skip the search step (only benchmark 3Di generation)
    uv run python benchmark_foldseek.py assemblies.txt \
        --foldseek-path foldseek/bin/foldseek \
        --prostt5-weights ~/prostT5-runner/prostt5_weights/prostt5/prostt5-f16.gguf \
        --gpu --thread-counts 1 4 8 --skip-search
"""

import argparse
import csv
import os
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

from Bio import SeqIO


def resolve_foldseek_bin(foldseek_path: str = None) -> str:
    if foldseek_path:
        if os.path.isdir(foldseek_path):
            bin_path = os.path.join(foldseek_path, "bin", "foldseek")
            if not os.path.exists(bin_path):
                bin_path = os.path.join(foldseek_path, "foldseek")
        else:
            bin_path = foldseek_path
        if not os.path.exists(bin_path):
            raise FileNotFoundError(f"foldseek binary not found at {bin_path}")
        return bin_path
    return "foldseek"


def measure_genome(
    fasta_path: str,
    foldseek_bin: str,
    prostt5_weights: str,
    threads: int,
    use_gpu: bool,
    skip_search: bool,
    work_dir: str,
) -> dict:
    """
    Run the foldseek pipeline on a single FASTA file and measure timing
    for each stage separately.

    Returns a dict with timing and size information.
    """
    db_prefix = os.path.join(work_dir, "benchDB")
    tmp_dir = os.path.join(work_dir, "tmp")
    result_prefix = os.path.join(work_dir, "result")
    result_file = os.path.join(work_dir, "result.tsv")

    # --- createdb (includes ProstT5 3Di prediction) ---
    cmd_create = [
        foldseek_bin, "createdb",
        fasta_path, db_prefix,
        "--prostt5-model", prostt5_weights,
        "--threads", str(threads),
    ]
    if use_gpu:
        cmd_create.extend(["--gpu", "1"])

    t0 = time.perf_counter()
    subprocess.run(cmd_create, capture_output=True, text=True, check=True)
    t_createdb = time.perf_counter() - t0

    # --- search (all-vs-all) ---
    t_search = 0.0
    if not skip_search:
        cmd_search = [
            foldseek_bin, "search",
            db_prefix, db_prefix,
            result_prefix, tmp_dir,
            "--threads", str(threads),
        ]
        t0 = time.perf_counter()
        subprocess.run(cmd_search, capture_output=True, text=True, check=True)
        t_search = time.perf_counter() - t0

        # convertalis (usually fast, include in search time)
        cmd_convert = [
            foldseek_bin, "convertalis",
            db_prefix, db_prefix,
            result_prefix, result_file,
        ]
        t0_conv = time.perf_counter()
        subprocess.run(cmd_convert, capture_output=True, text=True, check=True)
        t_search += time.perf_counter() - t0_conv

    return {
        "time_createdb_s": round(t_createdb, 3),
        "time_search_s": round(t_search, 3),
        "time_total_s": round(t_createdb + t_search, 3),
    }


def get_fasta_stats(fasta_path: str) -> tuple[int, int]:
    """Return (num_proteins, total_residues) for a FASTA file."""
    num_proteins = 0
    total_residues = 0
    for record in SeqIO.parse(fasta_path, "fasta"):
        num_proteins += 1
        total_residues += len(record.seq)
    return num_proteins, total_residues


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark foldseek 3Di prediction + search per genome across thread counts.",
    )
    parser.add_argument(
        "accession_file",
        help="Text file with one NCBI genome assembly accession per line",
    )
    parser.add_argument(
        "--foldseek-path", required=True,
        help="Path to foldseek binary",
    )
    parser.add_argument(
        "--prostt5-weights", required=True,
        help="Path to ProstT5 weights (e.g. prostt5-f16.gguf)",
    )
    parser.add_argument(
        "--gpu", action="store_true",
        help="Use GPU for ProstT5 inference",
    )
    parser.add_argument(
        "--thread-counts", nargs="+", type=int, default=[1, 2, 4, 8],
        help="Thread counts to benchmark (default: 1 2 4 8)",
    )
    parser.add_argument(
        "--reps", type=int, default=1,
        help="Number of repetitions per (genome, thread_count) combination (default: 1)",
    )
    parser.add_argument(
        "--fasta-dir", default=None,
        help=(
            "Directory containing pre-downloaded .faa files named {accession}.faa. "
            "If not provided, looks in batch_output/fastas/ then tries downloading."
        ),
    )
    parser.add_argument(
        "--skip-search", action="store_true",
        help="Only benchmark 3Di generation (createdb), skip the search step",
    )
    parser.add_argument(
        "--output", "-o", default="benchmark_results.csv",
        help="Output CSV file (default: benchmark_results.csv)",
    )
    args = parser.parse_args()

    foldseek_bin = resolve_foldseek_bin(args.foldseek_path)

    # Read accessions
    accessions = []
    with open(args.accession_file) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                accessions.append(line)

    # Locate FASTA files
    fasta_dirs = [
        args.fasta_dir,
        "batch_output/fastas",
        "results_fs_only_full_fs_compare/fastas",
    ]
    fasta_paths: dict[str, str] = {}
    for acc in accessions:
        for d in fasta_dirs:
            if d is None:
                continue
            candidate = os.path.join(d, f"{acc}.faa")
            if os.path.exists(candidate):
                fasta_paths[acc] = candidate
                break
        if acc not in fasta_paths:
            print(f"  [WARN] No .faa file found for {acc}, skipping")

    if not fasta_paths:
        print("No FASTA files found. Download them first with batch_3di_foldseek.py.")
        return

    # Compute genome sizes
    print("Measuring genome sizes...")
    genome_stats: dict[str, tuple[int, int]] = {}
    for acc, path in fasta_paths.items():
        n_prot, n_res = get_fasta_stats(path)
        genome_stats[acc] = (n_prot, n_res)
        print(f"  {acc}: {n_prot:,} proteins, {n_res:,} residues")

    # Sort by total residues (smallest first) for readable output
    sorted_accessions = sorted(fasta_paths.keys(), key=lambda a: genome_stats[a][1])

    # Run benchmarks
    results = []
    total_runs = len(sorted_accessions) * len(args.thread_counts) * args.reps
    run_num = 0

    for acc in sorted_accessions:
        fasta = fasta_paths[acc]
        n_prot, n_res = genome_stats[acc]

        for threads in args.thread_counts:
            for rep in range(1, args.reps + 1):
                run_num += 1
                label = f"[{run_num}/{total_runs}] {acc} threads={threads}"
                if args.reps > 1:
                    label += f" rep={rep}"
                print(f"\n{label} ({n_prot:,} proteins, {n_res:,} residues)")

                # Use a fresh temp directory each run to avoid caching effects
                with tempfile.TemporaryDirectory(prefix="bench_fs_") as work_dir:
                    try:
                        timing = measure_genome(
                            fasta, foldseek_bin, args.prostt5_weights,
                            threads, args.gpu, args.skip_search, work_dir,
                        )
                    except subprocess.CalledProcessError as e:
                        print(f"  FAILED: {e.stderr[:300] if e.stderr else e}")
                        continue

                row = {
                    "accession": acc,
                    "num_proteins": n_prot,
                    "total_residues": n_res,
                    "threads": threads,
                    "rep": rep,
                    **timing,
                }
                results.append(row)
                print(
                    f"  createdb: {timing['time_createdb_s']:.1f}s  "
                    f"search: {timing['time_search_s']:.1f}s  "
                    f"total: {timing['time_total_s']:.1f}s"
                )

    # Write CSV
    if results:
        fieldnames = [
            "accession", "num_proteins", "total_residues",
            "threads", "rep",
            "time_createdb_s", "time_search_s", "time_total_s",
        ]
        with open(args.output, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)
        print(f"\nResults written to {args.output}")

        # Print summary table
        print(f"\n{'='*80}")
        print(f"{'Accession':<20} {'Proteins':>10} {'Residues':>12} {'Threads':>8} {'CreateDB':>10} {'Search':>10} {'Total':>10}")
        print(f"{'-'*80}")
        for r in results:
            print(
                f"{r['accession']:<20} {r['num_proteins']:>10,} {r['total_residues']:>12,} "
                f"{r['threads']:>8} {r['time_createdb_s']:>9.1f}s {r['time_search_s']:>9.1f}s {r['time_total_s']:>9.1f}s"
            )
    else:
        print("\nNo results collected.")


if __name__ == "__main__":
    main()
