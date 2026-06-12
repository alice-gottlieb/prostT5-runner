#!/usr/bin/env python3
"""Build a target_id -> genome map by scanning per-genome FASTA (.faa) files.

foldseek_topn_pfam.py records each hit's target as a bare protein accession
(e.g. WP_000012948.1), but the matrix step needs to know which genome(s) that
protein belongs to. NCBI "WP_" proteins are multispecies, so one accession can
appear in many genomes; we therefore emit one row per (target_id, genome) pair
and let the downstream step count a hit toward every genome it appears in.

Input is a directory of FASTA files, one per genome, named <genome>.faa (e.g.
GCF_000005845.2.faa). Each header's first token after ">" is the protein
accession; the genome is the file's stem.

Pass --wanted-ids to restrict the map to a set of accessions (one per line) --
useful for a quick test against just the proteins a given run actually hit.
Without it, every protein in every file is mapped.

Output is a TSV with a header `target_id<TAB>genome`, matching the format
foldseek_topn_pfam.py's --target-genome-map expects.

Examples:
  # Map only the proteins this CPU test hit.
  uv run python build_target_genome_map.py \
      $SCRATCH/ncbi_genomes target_genome_map.tsv \
      --wanted-ids wanted_target_ids.txt

  # Map every protein in every genome (the full build).
  uv run python build_target_genome_map.py $SCRATCH/ncbi_genomes full_map.tsv
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime
from pathlib import Path

import polars as pl

# A byte that never appears in FASTA text, so scan_csv keeps each whole line in
# a single column instead of splitting it on commas/whitespace.
LINE_SEPARATOR = "\x01"


def log(message: str) -> None:
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}", flush=True)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("genomes_dir", type=Path,
                   help="Directory of per-genome FASTA files named <genome>.faa.")
    p.add_argument("output_tsv", type=Path,
                   help="Output TSV (target_id<TAB>genome, one row per pair).")
    p.add_argument("--wanted-ids", type=Path, default=None,
                   help="Optional file of protein accessions (one per line) to "
                        "keep; without it, all proteins are mapped.")
    p.add_argument("--glob", default="*.faa",
                   help="Filename glob for genome FASTAs (default: *.faa).")
    return p.parse_args()


def load_wanted_ids(path: Path | None) -> set[str] | None:
    if path is None:
        log("No --wanted-ids given; mapping every protein in every genome.")
        return None
    wanted = {line.strip() for line in path.read_text().splitlines() if line.strip()}
    log(f"Loaded {len(wanted):,} wanted target accessions from {path}.")
    return wanted


def build_map(genomes_dir: Path, wanted: set[str] | None, file_glob: str,
              out_tsv: Path) -> None:
    genome_files = sorted(genomes_dir.glob(file_glob))
    if not genome_files:
        sys.exit(f"ERROR: no files matching {file_glob} in {genomes_dir}")
    log(f"Scanning {len(genome_files):,} genome files in {genomes_dir}.")

    # Read every file at once as lines, tagged with the file they came from, then
    # keep only header lines and pull out the genome (file stem) and accession
    # (first non-whitespace token after ">"). One row per (target_id, genome).
    lines = pl.scan_csv(
        [str(f) for f in genome_files],
        has_header=False,
        separator=LINE_SEPARATOR,
        quote_char=None,
        new_columns=["line"],
        infer_schema=False,
        include_file_paths="path",
    )
    pairs = (
        lines
        .filter(pl.col("line").str.starts_with(">"))
        .select(
            pl.col("line").str.extract(r"^>(\S+)").alias("target_id"),
            pl.col("path").str.split("/").list.last()
                          .str.replace(r"\.[^.]+$", "").alias("genome"),
        )
    )
    if wanted is not None:
        pairs = pairs.filter(pl.col("target_id").is_in(list(wanted)))

    df = pairs.collect()
    df.write_csv(out_tsv, separator="\t")

    log(f"Wrote {df.height:,} (target, genome) pairs to {out_tsv}.")
    if wanted is not None:
        seen = set(df.get_column("target_id").unique().to_list())
        missing = wanted - seen
        log(f"Matched {len(seen):,}/{len(wanted):,} wanted accessions; "
            f"{len(missing):,} not found in any genome.")
        if missing:
            sample = ", ".join(sorted(missing)[:10])
            log(f"Examples of unmatched accessions: {sample}")


def main() -> int:
    args = parse_args()
    wanted = load_wanted_ids(args.wanted_ids)
    build_map(args.genomes_dir, wanted, args.glob, args.output_tsv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
