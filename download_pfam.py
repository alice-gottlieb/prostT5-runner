"""
Download Pfam family alignments from EBI InterPro and extract per-family
FASTA files suitable for foldseek 3Di prediction.

Downloads the bulk Stockholm-format alignment file (Pfam-A.seed.gz or
Pfam-A.full.gz) from the EBI FTP, parses each family block, and writes
one FASTA file per Pfam family with gap characters removed.

Output structure:
    <output_dir>/
        PF00001/PF00001.fasta
        PF00002/PF00002.fasta
        ...

Requires: requests, biopython (optional, only if --format=stockholm-bio)

Example usages:

    # Download seed alignments (smaller, faster)
    uv run python download_pfam.py --output-dir pfam_fastas

    # Download full alignments (larger, more sequences per family)
    uv run python download_pfam.py --output-dir pfam_fastas --alignment-type full

    # Resume interrupted download (skips already-extracted families)
    uv run python download_pfam.py --output-dir pfam_fastas --resume

    # Only download specific families
    uv run python download_pfam.py --output-dir pfam_fastas --families PF00001 PF00002 PF00076
"""

import argparse
import gzip
import os
import re
import sys
import time
from pathlib import Path
from urllib.request import urlretrieve, Request, urlopen
from urllib.error import URLError, HTTPError


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PFAM_FTP_BASE = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release"
PFAM_SEED_URL = f"{PFAM_FTP_BASE}/Pfam-A.seed.gz"
PFAM_FULL_URL = f"{PFAM_FTP_BASE}/Pfam-A.full.gz"


# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------

def download_with_retry(url: str, output_path: str, max_retries: int = 4) -> str:
    """Download a file with exponential backoff retry on failure."""
    if os.path.exists(output_path):
        print(f"  File already exists: {output_path}, skipping download")
        return output_path

    partial_path = output_path + ".partial"

    for attempt in range(max_retries):
        try:
            print(f"  Downloading {url}")
            print(f"  -> {output_path}")

            def _progress(block_num, block_size, total_size):
                downloaded = block_num * block_size
                if total_size > 0:
                    pct = min(100, downloaded * 100 / total_size)
                    mb = downloaded / (1024 * 1024)
                    total_mb = total_size / (1024 * 1024)
                    print(
                        f"\r  {mb:.1f} / {total_mb:.1f} MB ({pct:.1f}%)",
                        end="", flush=True,
                    )
                else:
                    mb = downloaded / (1024 * 1024)
                    print(f"\r  {mb:.1f} MB downloaded", end="", flush=True)

            urlretrieve(url, partial_path, reporthook=_progress)
            print()  # newline after progress
            os.rename(partial_path, output_path)
            print(f"  Download complete: {output_path}")
            return output_path

        except (URLError, HTTPError, OSError) as e:
            if os.path.exists(partial_path):
                os.remove(partial_path)
            if attempt < max_retries - 1:
                wait = 2 ** (attempt + 1)
                print(f"\n  Download failed: {e}. Retrying in {wait}s...")
                time.sleep(wait)
            else:
                print(f"\n  Download failed after {max_retries} attempts: {e}")
                raise


# ---------------------------------------------------------------------------
# Stockholm parser
# ---------------------------------------------------------------------------

def parse_stockholm_stream(fileobj, families_filter=None, resume_existing=None):
    """
    Parse a multi-Stockholm format file and yield per-family data.

    Each Stockholm alignment block starts with '# STOCKHOLM 1.0' and ends
    with '//'. Within each block:
      - '#=GF AC PFxxxxx.yy' gives the accession
      - '#=GF ID name' gives the family identifier
      - Sequence lines are: 'seqname/start-end  aligned_sequence'

    Yields:
        (accession, family_id, sequences)
        where sequences is a list of (seq_id, ungapped_sequence) tuples.
    """
    in_block = False
    accession = None
    family_id = None
    sequences = {}  # ordered dict of seqname -> aligned_seq (accumulated)

    for line in fileobj:
        if isinstance(line, bytes):
            line = line.decode("utf-8", errors="replace")
        line = line.rstrip("\n\r")

        if line.startswith("# STOCKHOLM 1.0"):
            in_block = True
            accession = None
            family_id = None
            sequences = {}
            continue

        if line == "//":
            if in_block and accession:
                # Strip version from accession (PF00001.23 -> PF00001)
                acc_base = accession.split(".")[0]

                # Skip if filtering and not in filter set
                if families_filter and acc_base not in families_filter:
                    in_block = False
                    continue

                # Skip if resuming and already exists
                if resume_existing and acc_base in resume_existing:
                    in_block = False
                    continue

                # Remove gap characters from sequences
                ungapped = []
                for seq_id, aligned_seq in sequences.items():
                    clean = aligned_seq.replace(".", "").replace("-", "")
                    if clean:  # skip empty sequences
                        ungapped.append((seq_id, clean))

                if ungapped:
                    yield (acc_base, family_id or acc_base, ungapped)

            in_block = False
            continue

        if not in_block:
            continue

        # Metadata lines
        if line.startswith("#=GF AC"):
            accession = line.split()[-1].strip()
        elif line.startswith("#=GF ID"):
            family_id = line.split(None, 2)[-1].strip()
        elif line.startswith("#") or not line.strip():
            continue
        else:
            # Sequence line: seqname  aligned_sequence
            parts = line.split()
            if len(parts) >= 2:
                seq_id = parts[0]
                aligned_seq = parts[-1]
                if seq_id in sequences:
                    sequences[seq_id] += aligned_seq
                else:
                    sequences[seq_id] = aligned_seq


def write_family_fasta(output_dir: str, accession: str, sequences: list) -> str:
    """Write a single family's sequences to a FASTA file."""
    family_dir = os.path.join(output_dir, accession)
    os.makedirs(family_dir, exist_ok=True)
    fasta_path = os.path.join(family_dir, f"{accession}.fasta")

    with open(fasta_path, "w") as f:
        for seq_id, seq in sequences:
            f.write(f">{seq_id}\n{seq}\n")

    return fasta_path


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def download_and_extract_pfam(
    output_dir: str,
    alignment_type: str = "seed",
    families_filter: set = None,
    resume: bool = False,
    download_dir: str = None,
):
    """
    Download the Pfam bulk alignment file and extract per-family FASTAs.

    Args:
        output_dir: Directory to write per-family FASTA files.
        alignment_type: 'seed' or 'full'.
        families_filter: Optional set of accessions (e.g. {'PF00001'}) to extract.
        resume: If True, skip families whose output directories already exist.
        download_dir: Directory for the raw downloaded .gz file. Defaults to output_dir.
    """
    url = PFAM_SEED_URL if alignment_type == "seed" else PFAM_FULL_URL
    filename = f"Pfam-A.{alignment_type}.gz"

    dl_dir = download_dir or output_dir
    os.makedirs(dl_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)

    gz_path = os.path.join(dl_dir, filename)

    # Download
    print(f"\n--- Downloading Pfam-A {alignment_type} alignments ---")
    download_with_retry(url, gz_path)

    # Determine which families to skip if resuming
    resume_existing = None
    if resume:
        resume_existing = set()
        if os.path.isdir(output_dir):
            for entry in os.listdir(output_dir):
                fasta_path = os.path.join(output_dir, entry, f"{entry}.fasta")
                if os.path.isfile(fasta_path):
                    resume_existing.add(entry)
        if resume_existing:
            print(f"  Resuming: skipping {len(resume_existing)} already-extracted families")

    # Parse and extract
    print(f"\n--- Parsing Stockholm file and extracting per-family FASTAs ---")
    family_count = 0
    seq_count = 0

    with gzip.open(gz_path, "rt", encoding="utf-8", errors="replace") as f:
        for accession, family_id, sequences in parse_stockholm_stream(
            f, families_filter=families_filter, resume_existing=resume_existing
        ):
            fasta_path = write_family_fasta(output_dir, accession, sequences)
            family_count += 1
            seq_count += len(sequences)

            if family_count % 500 == 0:
                print(f"  Extracted {family_count} families ({seq_count} sequences)...")

    print(f"\n--- Done ---")
    print(f"  Extracted {family_count} families, {seq_count} total sequences")
    print(f"  Output: {output_dir}/[PFxxxxx]/[PFxxxxx].fasta")

    return family_count


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Download Pfam family alignments from EBI and extract "
            "per-family FASTA files for foldseek 3Di prediction."
        )
    )
    parser.add_argument(
        "--output-dir", "-o",
        default="pfam_fastas",
        help="Directory for per-family FASTA output (default: pfam_fastas)",
    )
    parser.add_argument(
        "--download-dir",
        default=None,
        help="Directory for raw downloaded .gz file (default: same as output-dir)",
    )
    parser.add_argument(
        "--alignment-type",
        choices=["seed", "full"],
        default="seed",
        help="Which Pfam alignment to download: 'seed' (curated, smaller) "
             "or 'full' (all members, larger). Default: seed",
    )
    parser.add_argument(
        "--families",
        nargs="+",
        default=None,
        help="Only extract these Pfam accessions (e.g. PF00001 PF00002). "
             "Default: extract all.",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Skip families whose FASTA files already exist in the output directory.",
    )

    args = parser.parse_args()

    families_filter = None
    if args.families:
        families_filter = set(args.families)
        print(f"Filtering to {len(families_filter)} families: {families_filter}")

    download_and_extract_pfam(
        output_dir=args.output_dir,
        alignment_type=args.alignment_type,
        families_filter=families_filter,
        resume=args.resume,
        download_dir=args.download_dir,
    )


if __name__ == "__main__":
    main()
