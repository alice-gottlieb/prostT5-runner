#!/usr/bin/env python3
"""Check completion progress of 3DI chunks against completed accessions."""

import os
import subprocess
import sys
from pathlib import Path


def load_accessions_from_file(filepath):
    """Load accessions from a file (one per line, first whitespace-delimited field)."""
    accessions = set()
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line:
                accessions.add(line.split()[0])
    return accessions


def chunk_sort_key(filepath):
    """Extract numeric part from chunk filename for natural sorting."""
    import re
    name = filepath.name if isinstance(filepath, Path) else filepath
    numbers = re.findall(r'\d+', name)
    return (int(numbers[-1]) if numbers else 0, name)


def load_chunk_dir(chunk_dir):
    """Load all accessions from chunk files in a directory, sorted numerically."""
    accessions = {}
    files = [f for f in Path(chunk_dir).iterdir() if f.is_file()]
    for chunk_file in sorted(files, key=chunk_sort_key):
        with open(chunk_file) as f:
            accs = {line.strip() for line in f if line.strip()}
        accessions[chunk_file.name] = accs
    return accessions


def summarize_chunks(chunks, completed):
    """Group consecutive chunks by completion status and return summary lines."""
    chunk_names = list(chunks.keys())
    if not chunk_names:
        return []

    # Compute per-chunk percentage
    chunk_pcts = []
    for name in chunk_names:
        accs = chunks[name]
        total = len(accs)
        done = len(accs & completed)
        pct = (done / total * 100) if total else 0
        chunk_pcts.append((name, done, total, pct))

    def pct_bucket(pct):
        if pct == 100:
            return "100%"
        elif pct == 0:
            return "0%"
        else:
            return f"{int(pct // 10) * 10}-{int(pct // 10) * 10 + 10}%"

    # Group consecutive chunks with the same bucket
    lines = []
    i = 0
    while i < len(chunk_pcts):
        bucket = pct_bucket(chunk_pcts[i][3])
        j = i + 1
        group_done = chunk_pcts[i][1]
        group_total = chunk_pcts[i][2]
        while j < len(chunk_pcts) and pct_bucket(chunk_pcts[j][3]) == bucket:
            group_done += chunk_pcts[j][1]
            group_total += chunk_pcts[j][2]
            j += 1

        count = j - i
        if count == 1:
            name, done, total, pct = chunk_pcts[i]
            lines.append(f"  {name:<40} {done:>6}/{total:<6}  ({pct:.1f}%)")
        else:
            first = chunk_pcts[i][0]
            last = chunk_pcts[j - 1][0]
            group_pct = (group_done / group_total * 100) if group_total else 0
            lines.append(
                f"  {first} .. {last}  ({count} chunks)".ljust(42)
                + f"{group_done:>6}/{group_total:<6}  ({group_pct:.1f}%)"
            )
        i = j

    return lines


def get_chunk_pcts(chunks, completed):
    """Return list of per-chunk completion percentages."""
    pcts = []
    for name in chunks:
        accs = chunks[name]
        total = len(accs)
        done = len(accs & completed)
        pcts.append((done / total * 100) if total else 0)
    return pcts


def print_chunk_report(label, chunk_dir, completed):
    """Print a progress report for a chunk directory. Returns (label, pcts) or None."""
    if not Path(chunk_dir).is_dir():
        print(f"\n{'=' * 60}")
        print(f"  {label}: {chunk_dir}")
        print(f"  ERROR: directory does not exist")
        print(f"{'=' * 60}")
        return None

    chunks = load_chunk_dir(chunk_dir)
    all_in_chunks = set()
    for accs in chunks.values():
        all_in_chunks.update(accs)

    done = all_in_chunks & completed
    remaining = all_in_chunks - completed
    total = len(all_in_chunks)
    pct = (len(done) / total * 100) if total else 0

    n_complete = sum(1 for accs in chunks.values() if accs <= completed)
    n_partial = sum(1 for accs in chunks.values() if (accs & completed) and not (accs <= completed))
    n_empty = sum(1 for accs in chunks.values() if not (accs & completed))

    print(f"\n{'=' * 60}")
    print(f"  {label}")
    print(f"  {chunk_dir}")
    print(f"{'=' * 60}")
    print(f"  Total accessions:     {total:>8}")
    print(f"  Completed:            {len(done):>8}  ({pct:.1f}%)")
    print(f"  Remaining:            {len(remaining):>8}  ({100 - pct:.1f}%)")
    print(f"  Chunk files:          {len(chunks):>8}")
    print(f"    Fully done:         {n_complete:>8}")
    print(f"    Partially done:     {n_partial:>8}")
    print(f"    Not started:        {n_empty:>8}")
    print(f"{'─' * 60}")

    lines = summarize_chunks(chunks, completed)
    for line in lines:
        print(line)

    return (label, get_chunk_pcts(chunks, completed))


def plot_histograms(results, output_path):
    """Plot one histogram per folder showing chunk completion distribution."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(len(results), 1, figsize=(10, 5 * len(results)))
    if len(results) == 1:
        axes = [axes]

    bins = list(range(0, 110, 10))  # 0, 10, 20, ..., 100

    for ax, (label, pcts) in zip(axes, results):
        ax.hist(pcts, bins=bins, edgecolor="black", color="steelblue")
        ax.set_title(f"{label}  ({len(pcts)} chunks)")
        ax.set_xlabel("Chunk completion (%)")
        ax.set_ylabel("Number of chunks")
        ax.set_xticks(bins)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    print(f"\nHistogram saved to {output_path}")


def main():
    script_dir = Path(__file__).parent
    print("Running collect_accessions.py...")
    result = subprocess.run(
        ["uv", "run", "python", str(script_dir / "collect_accessions.py")],
        cwd=script_dir,
    )
    if result.returncode != 0:
        sys.exit("Error: collect_accessions.py failed")

    completed_file = script_dir / "all_accessions.txt"
    if not completed_file.is_file():
        sys.exit(f"Error: {completed_file} not found.")

    completed = load_accessions_from_file(completed_file)

    scratch = os.environ.get("SCRATCH")
    if not scratch:
        sys.exit("Error: $SCRATCH environment variable is not set")

    chunk_dirs = [
        ("all_3di_chunks_400_split", os.path.expanduser("~/prostT5-runner/all_3di_chunks_400_split")),
        ("new_chunks_20", os.path.join(scratch, "new_chunks_20")),
    ]

    print(f"\nCompleted accessions loaded: {len(completed)}")

    results = []
    for label, path in chunk_dirs:
        r = print_chunk_report(label, path, completed)
        if r:
            results.append(r)

    if results:
        plot_histograms(results, script_dir / "chunk_progress_histograms.png")

    print()


if __name__ == "__main__":
    main()
