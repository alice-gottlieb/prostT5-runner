"""
Build a counts matrix of Pfam domain 3Di hits across a species' genome.

Takes per-family Pfam 3Di FASTAs (from run_pfam_foldseek.sh) and a genome's
3Di FASTA (from batch_3di_foldseek.py), then scans each genome protein for
occurrences of each Pfam domain's 3Di signature using sliding-window matching.

Outputs:
  1. hits.tsv      — every hit with exact position, score, and sequences
  2. counts.tsv    — matrix of (genome_protein × pfam_domain) hit counts

Usage:
    # Basic usage
    uv run python pfam_3di_matrix.py \\
        --pfam-3di-dir $SCRATCH/pfam-3dis \\
        --genome-3di genome_3di.fasta \\
        --output-dir results

    # Stricter threshold (90% match required)
    uv run python pfam_3di_matrix.py \\
        --pfam-3di-dir $SCRATCH/pfam-3dis \\
        --genome-3di genome_3di.fasta \\
        --threshold 0.9

    # Use absolute mismatch count instead of fraction
    uv run python pfam_3di_matrix.py \\
        --pfam-3di-dir $SCRATCH/pfam-3dis \\
        --genome-3di genome_3di.fasta \\
        --max-mismatches 5

    # Run sanity tests
    uv run python pfam_3di_matrix.py --run-tests
"""

import argparse
import csv
import os
import sys
import tempfile
from collections import defaultdict
from pathlib import Path


# ---------------------------------------------------------------------------
# FASTA I/O
# ---------------------------------------------------------------------------

def parse_fasta(filepath: str) -> dict[str, str]:
    """Parse a FASTA file into {sequence_id: sequence}."""
    sequences = {}
    current_id = None
    current_seq: list[str] = []

    with open(filepath) as f:
        for line in f:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            elif line:
                current_seq.append(line)

    if current_id is not None:
        sequences[current_id] = "".join(current_seq)

    return sequences


# ---------------------------------------------------------------------------
# Representative selection
# ---------------------------------------------------------------------------

def select_representative(sequences: dict[str, str]) -> tuple[str, str]:
    """
    Select a representative sequence from a family by choosing the one
    closest to the median length. Ties broken by alphabetical ID.
    """
    if not sequences:
        raise ValueError("Cannot select representative from empty sequence set")

    items = sorted(sequences.items(), key=lambda x: len(x[1]))
    median_idx = len(items) // 2
    return items[median_idx]


# ---------------------------------------------------------------------------
# Sliding-window 3Di matching
# ---------------------------------------------------------------------------

def sliding_window_hits(
    query: str,
    target: str,
    threshold: float = 0.7,
    max_mismatches: int | None = None,
) -> list[dict]:
    """
    Slide query along target and return all positions where the match
    meets the threshold.

    Args:
        query: The Pfam domain 3Di sequence (the pattern to search for).
        target: The genome protein 3Di sequence (searched against).
        threshold: Minimum fraction of matching 3Di characters (0.0–1.0).
        max_mismatches: If set, overrides threshold with an absolute count.

    Returns:
        List of hit dicts: {position, match_fraction, mismatches, matched_region}
    """
    qlen = len(query)
    tlen = len(target)

    if qlen == 0 or qlen > tlen:
        return []

    query_lower = query.lower()
    target_lower = target.lower()

    hits = []
    n_windows = tlen - qlen + 1

    for i in range(n_windows):
        mismatches = 0
        mismatch_limit = max_mismatches if max_mismatches is not None else (qlen - int(qlen * threshold))

        for j in range(qlen):
            if query_lower[j] != target_lower[i + j]:
                mismatches += 1
                if mismatches > mismatch_limit:
                    break
        else:
            match_frac = (qlen - mismatches) / qlen
            if max_mismatches is not None or match_frac >= threshold:
                hits.append({
                    "position": i,
                    "match_fraction": round(match_frac, 4),
                    "mismatches": mismatches,
                    "matched_region": target[i : i + qlen],
                })

    return hits


# ---------------------------------------------------------------------------
# Load Pfam 3Di references
# ---------------------------------------------------------------------------

def load_pfam_representatives(
    pfam_3di_dir: str,
) -> dict[str, tuple[str, str]]:
    """
    Load one representative 3Di sequence per Pfam family from the
    directory structure produced by run_pfam_foldseek.sh.

    Expected layout:
        pfam_3di_dir/PFxxxxx/PFxxxxx_3di.fasta

    Returns:
        {pfam_accession: (representative_id, representative_3di_sequence)}
    """
    pfam_dir = Path(pfam_3di_dir)
    representatives = {}

    for family_dir in sorted(pfam_dir.iterdir()):
        if not family_dir.is_dir():
            continue

        pfam_acc = family_dir.name
        fasta_file = family_dir / f"{pfam_acc}_3di.fasta"

        if not fasta_file.exists():
            continue

        seqs = parse_fasta(str(fasta_file))
        if not seqs:
            continue

        rep_id, rep_seq = select_representative(seqs)
        representatives[pfam_acc] = (rep_id, rep_seq)

    return representatives


# ---------------------------------------------------------------------------
# Matrix construction
# ---------------------------------------------------------------------------

def scan_genome(
    pfam_representatives: dict[str, tuple[str, str]],
    genome_3di: dict[str, str],
    threshold: float = 0.7,
    max_mismatches: int | None = None,
) -> tuple[list[dict], dict[str, dict[str, int]]]:
    """
    Scan all genome proteins against all Pfam domain representatives.

    Returns:
        (all_hits, counts_matrix)

        all_hits: List of dicts with keys:
            protein_id, pfam_id, pfam_rep_id, position, match_fraction,
            mismatches, query_length, matched_region

        counts_matrix: {protein_id: {pfam_id: hit_count}}
    """
    all_hits = []
    counts: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))

    n_proteins = len(genome_3di)
    n_pfam = len(pfam_representatives)
    total_comparisons = n_proteins * n_pfam
    comparisons_done = 0
    last_pct = -1

    for protein_id, protein_seq in genome_3di.items():
        for pfam_acc, (rep_id, rep_seq) in pfam_representatives.items():
            comparisons_done += 1

            hits = sliding_window_hits(
                query=rep_seq,
                target=protein_seq,
                threshold=threshold,
                max_mismatches=max_mismatches,
            )

            if hits:
                counts[protein_id][pfam_acc] += len(hits)

                for h in hits:
                    all_hits.append({
                        "protein_id": protein_id,
                        "pfam_id": pfam_acc,
                        "pfam_rep_id": rep_id,
                        "position": h["position"],
                        "match_fraction": h["match_fraction"],
                        "mismatches": h["mismatches"],
                        "query_length": len(rep_seq),
                        "matched_region": h["matched_region"],
                    })

            pct = (comparisons_done * 100) // total_comparisons
            if pct >= last_pct + 5:
                last_pct = pct
                print(
                    f"  {comparisons_done}/{total_comparisons} comparisons "
                    f"({pct}%), {len(all_hits)} hits so far"
                )

    return all_hits, dict(counts)


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

HIT_COLUMNS = [
    "protein_id", "pfam_id", "pfam_rep_id", "position",
    "match_fraction", "mismatches", "query_length", "matched_region",
]


def write_hits_tsv(hits: list[dict], output_path: str):
    """Write detailed hits to a TSV file."""
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=HIT_COLUMNS, delimiter="\t")
        writer.writeheader()
        for hit in hits:
            writer.writerow(hit)
    print(f"  Wrote {len(hits)} hits to {output_path}")


def write_counts_matrix(
    counts: dict[str, dict[str, int]],
    all_pfam_ids: list[str],
    all_protein_ids: list[str],
    output_path: str,
):
    """Write the full protein × Pfam counts matrix to a TSV file."""
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["protein_id"] + all_pfam_ids)

        for protein_id in all_protein_ids:
            row = [protein_id]
            protein_counts = counts.get(protein_id, {})
            for pfam_id in all_pfam_ids:
                row.append(protein_counts.get(pfam_id, 0))
            writer.writerow(row)

    print(f"  Wrote {len(all_protein_ids)} × {len(all_pfam_ids)} matrix to {output_path}")


# ---------------------------------------------------------------------------
# Sanity tests
# ---------------------------------------------------------------------------

def _write_temp_fasta(seqs: dict[str, str], tmpdir: str, name: str) -> str:
    path = os.path.join(tmpdir, name)
    with open(path, "w") as f:
        for sid, seq in seqs.items():
            f.write(f">{sid}\n{seq}\n")
    return path


def _assert(condition, msg):
    if not condition:
        print(f"  FAIL: {msg}")
        sys.exit(1)
    print(f"  PASS: {msg}")


def run_sanity_tests():
    """Run all embedded sanity tests."""
    print("\n=== Sanity Tests ===\n")

    # --- Test 1: parse_fasta ---
    print("Test 1: parse_fasta")
    with tempfile.TemporaryDirectory() as tmpdir:
        path = _write_temp_fasta(
            {"seq1": "acdefgh", "seq2": "klmnpq"}, tmpdir, "test.fasta"
        )
        result = parse_fasta(path)
        _assert(len(result) == 2, "parsed 2 sequences")
        _assert(result["seq1"] == "acdefgh", "seq1 content correct")
        _assert(result["seq2"] == "klmnpq", "seq2 content correct")

    # Multi-line sequences
    print("\nTest 1b: parse_fasta multi-line")
    with tempfile.TemporaryDirectory() as tmpdir:
        path = os.path.join(tmpdir, "multi.fasta")
        with open(path, "w") as f:
            f.write(">s1\nacde\nfghi\n>s2\nklmn\n")
        result = parse_fasta(path)
        _assert(result["s1"] == "acdefghi", "multi-line concatenated")
        _assert(result["s2"] == "klmn", "single-line still works")

    # --- Test 2: select_representative ---
    print("\nTest 2: select_representative")
    seqs = {"short": "abc", "medium": "abcde", "long": "abcdefgh"}
    rep_id, rep_seq = select_representative(seqs)
    _assert(rep_id == "medium", f"median-length selected (got {rep_id})")
    _assert(rep_seq == "abcde", "correct sequence returned")

    seqs_single = {"only": "abcdef"}
    rep_id, _ = select_representative(seqs_single)
    _assert(rep_id == "only", "single sequence returns itself")

    # --- Test 3: sliding_window_hits (exact match) ---
    print("\nTest 3: sliding_window_hits — exact match")
    hits = sliding_window_hits("acdef", "xxacdefxx", threshold=1.0)
    _assert(len(hits) == 1, f"found 1 exact hit (got {len(hits)})")
    _assert(hits[0]["position"] == 2, f"at position 2 (got {hits[0]['position']})")
    _assert(hits[0]["match_fraction"] == 1.0, "100% match")
    _assert(hits[0]["mismatches"] == 0, "0 mismatches")

    # --- Test 4: sliding_window_hits (partial match with threshold) ---
    print("\nTest 4: sliding_window_hits — partial match")
    # query: "abcde", target has "aXcde" at pos 0 → 4/5 = 0.8 match
    hits = sliding_window_hits("abcde", "aXcdefffff", threshold=0.7)
    _assert(len(hits) >= 1, f"found hit at 70% threshold (got {len(hits)})")
    _assert(hits[0]["position"] == 0, "hit at position 0")
    _assert(hits[0]["match_fraction"] == 0.8, f"80% match (got {hits[0]['match_fraction']})")

    # Same query at stricter threshold
    hits_strict = sliding_window_hits("abcde", "aXcdefffff", threshold=0.9)
    _assert(len(hits_strict) == 0, "no hit at 90% threshold")

    # --- Test 5: sliding_window_hits (max_mismatches mode) ---
    print("\nTest 5: sliding_window_hits — max_mismatches")
    hits = sliding_window_hits("abcde", "aXcdefffff", max_mismatches=1)
    _assert(len(hits) >= 1, "found hit with max_mismatches=1")

    hits_zero = sliding_window_hits("abcde", "aXcdefffff", max_mismatches=0)
    _assert(len(hits_zero) == 0, "no hit with max_mismatches=0")

    # --- Test 6: sliding_window_hits (case insensitive) ---
    print("\nTest 6: sliding_window_hits — case insensitive")
    hits = sliding_window_hits("ACDEF", "xxacdefxx", threshold=1.0)
    _assert(len(hits) == 1, "case-insensitive match works")

    # --- Test 7: sliding_window_hits (query longer than target) ---
    print("\nTest 7: sliding_window_hits — query longer than target")
    hits = sliding_window_hits("abcdefghij", "abc", threshold=0.5)
    _assert(len(hits) == 0, "no hits when query > target")

    # --- Test 8: sliding_window_hits (multiple hits) ---
    print("\nTest 8: sliding_window_hits — multiple hits")
    hits = sliding_window_hits("abc", "abcXXabc", threshold=1.0)
    _assert(len(hits) == 2, f"found 2 exact hits (got {len(hits)})")
    _assert(hits[0]["position"] == 0, "first at 0")
    _assert(hits[1]["position"] == 5, f"second at 5 (got {hits[1]['position']})")

    # --- Test 9: sliding_window_hits (early termination) ---
    print("\nTest 9: early termination doesn't miss hits")
    # A long target where the hit is at the end
    target = "x" * 100 + "abcde"
    hits = sliding_window_hits("abcde", target, threshold=1.0)
    _assert(len(hits) == 1, "found hit at end of long target")
    _assert(hits[0]["position"] == 100, f"at position 100 (got {hits[0]['position']})")

    # --- Test 10: full pipeline integration ---
    print("\nTest 10: full pipeline integration")
    with tempfile.TemporaryDirectory() as tmpdir:
        # Set up fake Pfam 3Di directory
        pfam_dir = os.path.join(tmpdir, "pfam_3dis")
        for fam, seq in [("PF00001", "acdefghik"), ("PF00002", "lmnpqrstv")]:
            fam_dir = os.path.join(pfam_dir, fam)
            os.makedirs(fam_dir)
            _write_temp_fasta(
                {f"{fam}_rep": seq}, fam_dir, f"{fam}_3di.fasta"
            )

        # Set up fake genome 3Di
        genome_path = _write_temp_fasta(
            {
                "protein_A": "xxxacdefghikxxx",   # contains PF00001 at pos 3
                "protein_B": "lmnpqrstvxxxxxx",    # contains PF00002 at pos 0
                "protein_C": "xxxxxxxxxxxxxxxxx",  # no hits
            },
            tmpdir,
            "genome_3di.fasta",
        )

        # Load and scan
        reps = load_pfam_representatives(pfam_dir)
        _assert(len(reps) == 2, f"loaded 2 Pfam families (got {len(reps)})")

        genome = parse_fasta(genome_path)
        all_hits, counts = scan_genome(reps, genome, threshold=1.0)

        _assert(len(all_hits) == 2, f"found 2 total hits (got {len(all_hits)})")

        # protein_A should have PF00001
        pA_hits = [h for h in all_hits if h["protein_id"] == "protein_A"]
        _assert(len(pA_hits) == 1, "protein_A has 1 hit")
        _assert(pA_hits[0]["pfam_id"] == "PF00001", "protein_A hit is PF00001")
        _assert(pA_hits[0]["position"] == 3, f"at position 3 (got {pA_hits[0]['position']})")

        # protein_B should have PF00002
        pB_hits = [h for h in all_hits if h["protein_id"] == "protein_B"]
        _assert(len(pB_hits) == 1, "protein_B has 1 hit")
        _assert(pB_hits[0]["pfam_id"] == "PF00002", "protein_B hit is PF00002")
        _assert(pB_hits[0]["position"] == 0, f"at position 0 (got {pB_hits[0]['position']})")

        # protein_C should have no hits
        _assert("protein_C" not in counts, "protein_C has no hits in counts")

        # Counts matrix values
        _assert(counts["protein_A"]["PF00001"] == 1, "counts[A][PF00001] == 1")
        _assert(counts["protein_B"]["PF00002"] == 1, "counts[B][PF00002] == 1")

        # Write outputs and verify files
        out_dir = os.path.join(tmpdir, "output")
        os.makedirs(out_dir)

        hits_path = os.path.join(out_dir, "hits.tsv")
        write_hits_tsv(all_hits, hits_path)
        _assert(os.path.exists(hits_path), "hits.tsv created")

        with open(hits_path) as f:
            lines = f.readlines()
        _assert(len(lines) == 3, f"hits.tsv has header + 2 rows (got {len(lines)})")

        matrix_path = os.path.join(out_dir, "counts.tsv")
        all_pfam_ids = sorted(reps.keys())
        all_protein_ids = sorted(genome.keys())
        write_counts_matrix(counts, all_pfam_ids, all_protein_ids, matrix_path)
        _assert(os.path.exists(matrix_path), "counts.tsv created")

        with open(matrix_path) as f:
            lines = f.readlines()
        _assert(
            len(lines) == 4,
            f"counts.tsv has header + 3 protein rows (got {len(lines)})",
        )

        # Verify matrix content
        with open(matrix_path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)
        pC_row = [r for r in rows if r["protein_id"] == "protein_C"][0]
        _assert(
            pC_row["PF00001"] == "0" and pC_row["PF00002"] == "0",
            "protein_C row is all zeros",
        )

    # --- Test 11: threshold boundary ---
    print("\nTest 11: threshold boundary")
    # 4/5 = 0.8 exactly
    hits = sliding_window_hits("abcde", "abcXe", threshold=0.8)
    _assert(len(hits) == 1, "hit at exact threshold boundary (0.8)")
    hits = sliding_window_hits("abcde", "abcXe", threshold=0.81)
    _assert(len(hits) == 0, "no hit just above threshold (0.81)")

    # --- Test 12: overlapping hits ---
    print("\nTest 12: overlapping hits reported correctly")
    hits = sliding_window_hits("aaa", "aaaa", threshold=1.0)
    _assert(len(hits) == 2, f"2 overlapping hits in 'aaaa' (got {len(hits)})")
    _assert(hits[0]["position"] == 0, "first at 0")
    _assert(hits[1]["position"] == 1, "second at 1")

    # --- Test 13: empty inputs ---
    print("\nTest 13: empty/edge-case inputs")
    hits = sliding_window_hits("", "abcde", threshold=0.5)
    _assert(len(hits) == 0, "empty query yields no hits")
    hits = sliding_window_hits("abc", "", threshold=0.5)
    _assert(len(hits) == 0, "empty target yields no hits")

    print("\n=== All sanity tests passed ===\n")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Build a counts matrix of Pfam domain 3Di hits across a genome. "
            "Scans genome protein 3Di sequences for Pfam domain 3Di patterns "
            "using sliding-window matching."
        ),
    )
    parser.add_argument(
        "--pfam-3di-dir",
        required="--run-tests" not in sys.argv,
        help=(
            "Directory of per-family 3Di FASTAs from run_pfam_foldseek.sh. "
            "Expected structure: <dir>/PFxxxxx/PFxxxxx_3di.fasta"
        ),
    )
    parser.add_argument(
        "--genome-3di",
        required="--run-tests" not in sys.argv,
        help="FASTA file of 3Di sequences for genome proteins",
    )
    parser.add_argument(
        "--output-dir", "-o",
        default="pfam_matrix_output",
        help="Directory for output files (default: pfam_matrix_output)",
    )
    parser.add_argument(
        "--threshold", "-t",
        type=float,
        default=0.7,
        help="Minimum fraction of matching 3Di characters for a hit (default: 0.7)",
    )
    parser.add_argument(
        "--max-mismatches",
        type=int,
        default=None,
        help=(
            "Maximum allowed mismatches (absolute count). "
            "Overrides --threshold when set."
        ),
    )
    parser.add_argument(
        "--run-tests",
        action="store_true",
        help="Run embedded sanity tests and exit",
    )

    args = parser.parse_args()

    if args.run_tests:
        run_sanity_tests()
        return

    # Validate inputs
    if not os.path.isdir(args.pfam_3di_dir):
        print(f"ERROR: Pfam 3Di directory not found: {args.pfam_3di_dir}")
        sys.exit(1)
    if not os.path.isfile(args.genome_3di):
        print(f"ERROR: Genome 3Di FASTA not found: {args.genome_3di}")
        sys.exit(1)

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # 1. Load Pfam representatives
    print("\n--- Loading Pfam 3Di representatives ---")
    pfam_reps = load_pfam_representatives(args.pfam_3di_dir)
    print(f"  Loaded {len(pfam_reps)} Pfam families")
    if not pfam_reps:
        print("  No Pfam families found. Check --pfam-3di-dir structure.")
        sys.exit(1)

    # Report domain length stats
    lengths = [len(seq) for _, seq in pfam_reps.values()]
    print(f"  Domain lengths: min={min(lengths)}, median={sorted(lengths)[len(lengths)//2]}, max={max(lengths)}")

    # 2. Load genome 3Di sequences
    print("\n--- Loading genome 3Di sequences ---")
    genome_3di = parse_fasta(args.genome_3di)
    print(f"  Loaded {len(genome_3di)} proteins")
    if not genome_3di:
        print("  No sequences found in genome FASTA.")
        sys.exit(1)

    # 3. Scan
    mode = (
        f"max_mismatches={args.max_mismatches}"
        if args.max_mismatches is not None
        else f"threshold={args.threshold}"
    )
    print(f"\n--- Scanning genome × Pfam ({mode}) ---")
    all_hits, counts = scan_genome(
        pfam_reps,
        genome_3di,
        threshold=args.threshold,
        max_mismatches=args.max_mismatches,
    )

    # 4. Write outputs
    print("\n--- Writing outputs ---")
    all_pfam_ids = sorted(pfam_reps.keys())
    all_protein_ids = sorted(genome_3di.keys())

    hits_path = str(out / "hits.tsv")
    write_hits_tsv(all_hits, hits_path)

    matrix_path = str(out / "counts.tsv")
    write_counts_matrix(counts, all_pfam_ids, all_protein_ids, matrix_path)

    # 5. Summary
    proteins_with_hits = len(counts)
    pfam_with_hits = len({h["pfam_id"] for h in all_hits})
    print(f"\n--- Summary ---")
    print(f"  Total hits:               {len(all_hits)}")
    print(f"  Proteins with >= 1 hit:   {proteins_with_hits}/{len(genome_3di)}")
    print(f"  Pfam domains with hits:   {pfam_with_hits}/{len(pfam_reps)}")
    print(f"  Output:                   {out}/")


if __name__ == "__main__":
    main()
