"""
Analyze foldseek all-vs-all results and inspect 3Di codes.

Usage examples:

    # Overview stats for a results directory
    uv run python analyze_results.py results_fs_only_full_fs_compare overview

    # Show 3Di code for a specific protein
    uv run python analyze_results.py results_fs_only_full_fs_compare lookup WP_004995605.1

    # Find proteins similar to a target (sorted by e-value)
    uv run python analyze_results.py results_fs_only_full_fs_compare similar WP_004995605.1

    # Find similar proteins with e-value cutoff
    uv run python analyze_results.py results_fs_only_full_fs_compare similar WP_004995605.1 --evalue 1e-10

    # Export similar proteins to FASTA (AA sequences)
    uv run python analyze_results.py results_fs_only_full_fs_compare similar WP_004995605.1 --evalue 1e-10 --output-fasta hits.fasta
"""

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

import pandas as pd
from Bio import SeqIO


RESULT_COLUMNS = [
    "query", "target", "fident", "alnlen", "mismatch", "gapopen",
    "qstart", "qend", "tstart", "tend", "evalue", "bits",
]


def load_results(results_dir: Path) -> pd.DataFrame:
    tsv = results_dir / "foldseek_results" / "all_vs_all_results.tsv"
    df = pd.read_csv(tsv, sep="\t", header=None, names=RESULT_COLUMNS)
    return df


def load_3di_index(results_dir: Path) -> dict[str, str]:
    """Load 3Di FASTA into a {seq_id: 3di_sequence} dict."""
    fasta = results_dir / "foldseek_db" / "all_sequences_3di.fasta"
    return {r.id: str(r.seq) for r in SeqIO.parse(str(fasta), "fasta")}


def load_aa_index(results_dir: Path) -> dict[str, str]:
    """Load combined AA FASTA into a {seq_id: aa_sequence} dict."""
    fasta = results_dir / "combined_aa.fasta"
    return {r.id: str(r.seq) for r in SeqIO.parse(str(fasta), "fasta")}


def cmd_overview(results_dir: Path):
    """Print summary statistics about the database and search results."""
    meta_path = results_dir / "metadata.json"
    if meta_path.exists():
        meta = json.loads(meta_path.read_text())
        print("=== Metadata ===")
        print(f"  Assemblies requested:  {len(meta['assemblies_requested'])}")
        print(f"  Assemblies downloaded: {len(meta['assemblies_downloaded'])}")
        print(f"  Total sequences:       {meta['num_sequences']}")
        print()

    three_di = load_3di_index(results_dir)
    print(f"=== 3Di Database ===")
    print(f"  Sequences in 3Di FASTA: {len(three_di)}")
    lengths = [len(s) for s in three_di.values()]
    print(f"  Sequence length range:  {min(lengths)} - {max(lengths)}")
    print(f"  Mean length:            {sum(lengths)/len(lengths):.0f}")
    print()

    # Check for empty/trivial 3Di predictions
    empty = sum(1 for s in three_di.values() if len(s) < 5)
    all_same = sum(1 for s in three_di.values() if len(set(s)) <= 2)
    print(f"  Very short (<5 residues): {empty}")
    print(f"  Low complexity (<=2 unique chars): {all_same}")

    # 3Di alphabet usage
    all_chars = "".join(three_di.values())
    char_counts = Counter(all_chars)
    print(f"  3Di alphabet usage ({len(char_counts)} unique chars):")
    for ch, cnt in char_counts.most_common():
        print(f"    {ch}: {cnt:>10} ({100*cnt/len(all_chars):.1f}%)")
    print()

    print("=== Search Results ===")
    df = load_results(results_dir)
    print(f"  Total alignments:       {len(df):,}")
    # Exclude self-hits
    non_self = df[df["query"] != df["target"]]
    print(f"  Non-self alignments:    {len(non_self):,}")
    print(f"  Unique query proteins:  {df['query'].nunique():,}")
    print(f"  Unique target proteins: {df['target'].nunique():,}")
    print()

    print("  E-value distribution (non-self hits):")
    thresholds = [1e-50, 1e-30, 1e-10, 1e-5, 1e-3, 0.01, 0.1, 1.0]
    for t in thresholds:
        n = (non_self["evalue"] <= t).sum()
        print(f"    e <= {t:.0e}: {n:>10,} ({100*n/len(non_self):.1f}%)")
    print()

    print("  Sequence identity distribution (non-self hits):")
    for lo, hi in [(0.9, 1.0), (0.7, 0.9), (0.5, 0.7), (0.3, 0.5), (0.0, 0.3)]:
        n = ((non_self["fident"] >= lo) & (non_self["fident"] < hi)).sum()
        print(f"    [{lo:.1f}, {hi:.1f}): {n:>10,} ({100*n/len(non_self):.1f}%)")
    print()

    # Proteins with most hits
    hits_per_query = non_self.groupby("query").size().sort_values(ascending=False)
    print("  Top 10 proteins by number of hits:")
    for pid, count in hits_per_query.head(10).items():
        print(f"    {pid}: {count} hits")


def cmd_lookup(results_dir: Path, protein_id: str):
    """Show the 3Di code and AA sequence for a specific protein."""
    three_di = load_3di_index(results_dir)
    aa = load_aa_index(results_dir)

    if protein_id not in three_di:
        # Fuzzy match
        matches = [k for k in three_di if protein_id in k]
        if matches:
            print(f"Protein '{protein_id}' not found. Did you mean:")
            for m in matches[:10]:
                print(f"  {m}")
        else:
            print(f"Protein '{protein_id}' not found in 3Di database.")
            print(f"Database contains {len(three_di)} sequences.")
        return

    tdi = three_di[protein_id]
    aa_seq = aa.get(protein_id, "(not found in combined AA FASTA)")

    print(f"=== {protein_id} ===")
    print(f"  AA length:  {len(aa_seq) if isinstance(aa_seq, str) and aa_seq[0] != '(' else 'N/A'}")
    print(f"  3Di length: {len(tdi)}")
    print()
    print("AA sequence:")
    print(aa_seq)
    print()
    print("3Di sequence:")
    print(tdi)
    print()

    # 3Di composition
    char_counts = Counter(tdi)
    print("3Di composition:")
    for ch, cnt in char_counts.most_common():
        print(f"  {ch}: {cnt} ({100*cnt/len(tdi):.1f}%)")


def cmd_similar(results_dir: Path, protein_id: str, evalue: float, max_hits: int, output_fasta: str | None):
    """Find proteins similar to a target, sorted by e-value."""
    df = load_results(results_dir)

    hits = df[df["query"] == protein_id].copy()
    if hits.empty:
        # Check if it exists at all
        three_di = load_3di_index(results_dir)
        if protein_id not in three_di:
            matches = [k for k in three_di if protein_id in k]
            if matches:
                print(f"Protein '{protein_id}' not found. Did you mean:")
                for m in matches[:10]:
                    print(f"  {m}")
            else:
                print(f"Protein '{protein_id}' not found in database.")
            return
        print(f"No hits found for {protein_id} (it exists but has no alignments).")
        return

    # Filter out self-hit and apply e-value threshold
    hits = hits[hits["target"] != protein_id]
    hits = hits[hits["evalue"] <= evalue]
    hits = hits.sort_values("evalue")
    hits = hits.head(max_hits)

    if hits.empty:
        print(f"No hits for {protein_id} at e-value <= {evalue:.0e}")
        return

    print(f"=== {len(hits)} hits for {protein_id} (e-value <= {evalue:.0e}) ===\n")
    print(f"{'target':<25} {'fident':>7} {'alnlen':>7} {'evalue':>12} {'bits':>7}")
    print("-" * 62)
    for _, row in hits.iterrows():
        print(f"{row['target']:<25} {row['fident']:>7.3f} {row['alnlen']:>7} {row['evalue']:>12.2e} {row['bits']:>7}")

    if output_fasta:
        aa = load_aa_index(results_dir)
        with open(output_fasta, "w") as f:
            for _, row in hits.iterrows():
                tid = row["target"]
                if tid in aa:
                    f.write(f">{tid} evalue={row['evalue']:.2e} fident={row['fident']:.3f}\n")
                    f.write(aa[tid] + "\n")
        print(f"\nWrote {len(hits)} AA sequences to {output_fasta}")


def main():
    parser = argparse.ArgumentParser(description="Analyze foldseek results")
    parser.add_argument("results_dir", type=Path, help="Path to results directory")
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("overview", help="Summary statistics")

    p_lookup = sub.add_parser("lookup", help="Show 3Di/AA for a protein")
    p_lookup.add_argument("protein_id", help="Protein accession (e.g. WP_004995605.1)")

    p_sim = sub.add_parser("similar", help="Find similar proteins")
    p_sim.add_argument("protein_id", help="Query protein accession")
    p_sim.add_argument("--evalue", type=float, default=0.01, help="E-value cutoff (default: 0.01)")
    p_sim.add_argument("--max-hits", type=int, default=50, help="Max hits to show (default: 50)")
    p_sim.add_argument("--output-fasta", help="Write hit AA sequences to this FASTA file")

    args = parser.parse_args()

    if args.command == "overview":
        cmd_overview(args.results_dir)
    elif args.command == "lookup":
        cmd_lookup(args.results_dir, args.protein_id)
    elif args.command == "similar":
        cmd_similar(args.results_dir, args.protein_id, args.evalue, args.max_hits, args.output_fasta)


if __name__ == "__main__":
    main()
