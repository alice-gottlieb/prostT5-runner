"""
Batch download whole-species protein FASTA (.faa) files from NCBI Datasets,
convert all sequences to 3Di codes using foldseek's integrated ProstT5, then
run an all-vs-all foldseek search across all sequences.

This is an alternative to batch_3di_search.py that uses foldseek for both
3Di prediction and search, rather than independently loading ProstT5 from
HuggingFace via Python.

Input:  a text file with one NCBI genome assembly accession (GCF_* or GCA_*)
        per line
Output: per-assembly .faa files, a combined 3Di FASTA, and foldseek results

Requires: foldseek (v9+ with ProstT5 support)

Example usages:

    # Basic usage (auto-downloads ProstT5 weights, uses CPU)
    uv run python batch_3di_foldseek.py assemblies.txt

    # Custom output directory
    uv run python batch_3di_foldseek.py assemblies.txt -o my_results

    # With GPU acceleration
    uv run python batch_3di_foldseek.py assemblies.txt --gpu

    # Multi-threaded with GPU
    uv run python batch_3di_foldseek.py assemblies.txt --gpu --threads 8

    # With pre-downloaded ProstT5 weights (skips download on subsequent runs)
    uv run python batch_3di_foldseek.py assemblies.txt --prostt5-weights /path/to/prostt5

    # With a specific foldseek binary
    uv run python batch_3di_foldseek.py assemblies.txt --foldseek-path /opt/foldseek/bin/foldseek

    # With NCBI API key (faster downloads, 10 req/s vs 3 req/s)
    uv run python batch_3di_foldseek.py assemblies.txt --ncbi-api-key YOUR_KEY

    # Only generate 3Di codes, skip the all-vs-all search
    uv run python batch_3di_foldseek.py assemblies.txt --skip-foldseek

    # Pass extra arguments to foldseek search (e.g. e-value cutoff)
    uv run python batch_3di_foldseek.py assemblies.txt --foldseek-args -e 0.001

    # Full example with multiple options
    uv run python batch_3di_foldseek.py assemblies.txt \\
        -o results --gpu --threads 8 \\
        --ncbi-api-key YOUR_KEY \\
        --foldseek-args -e 0.001 --max-seqs 1000
"""

import io
import time
import json
import os
import subprocess
import argparse
import zipfile
import requests
from pathlib import Path
from Bio import SeqIO


# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------

def read_accessions(accession_file: str) -> list[str]:
    """Read genome assembly accessions from a text file, one per line.
    Ignores blank lines and comments (#)."""
    accessions = []
    with open(accession_file) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                accessions.append(line)
    print(f"Loaded {len(accessions)} assembly accessions from {accession_file}")
    return accessions


def download_species_faa(
    assembly_accession: str,
    output_path: str,
    ncbi_api_key: str = None,
) -> str | None:
    """
    Download the full protein FASTA (.faa) for a genome assembly using the
    NCBI Datasets API v2.

    The API returns a ZIP whose relevant entry is:
      ncbi_dataset/data/<accession>/protein.faa

    Args:
        assembly_accession: A RefSeq (GCF_*) or GenBank (GCA_*) assembly accession.
        output_path: Where to write the extracted .faa file.
        ncbi_api_key: Optional NCBI API key (increases rate limit to 10 req/s).

    Returns:
        output_path on success, None on failure.
    """
    url = (
        "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"
        f"{assembly_accession}/download"
    )
    params = {"include_annotation_type": "PROT_FASTA"}
    headers = {"Accept": "application/zip"}
    if ncbi_api_key:
        headers["api-key"] = ncbi_api_key

    try:
        resp = requests.get(url, params=params, headers=headers, timeout=120, stream=True)
        resp.raise_for_status()

        # Load the ZIP from memory and extract the protein FASTA
        zip_bytes = io.BytesIO(resp.content)
        with zipfile.ZipFile(zip_bytes) as zf:
            # Find the protein.faa entry (may be nested under the accession dir)
            faa_entries = [
                n for n in zf.namelist()
                if n.endswith("protein.faa")
            ]
            if not faa_entries:
                print(f"  [WARN] {assembly_accession}: no protein.faa found in downloaded zip")
                return None

            with zf.open(faa_entries[0]) as src, open(output_path, "wb") as dst:
                dst.write(src.read())

        return output_path

    except requests.RequestException as e:
        print(f"  [WARN] {assembly_accession}: download failed — {e}")
        return None
    except zipfile.BadZipFile as e:
        print(f"  [WARN] {assembly_accession}: bad zip response — {e}")
        return None


def batch_download(
    accessions: list[str],
    fasta_dir: str,
    ncbi_api_key: str = None,
    delay: float = 0.34,
) -> dict[str, str]:
    """
    Download .faa files for all genome assembly accessions.

    Respects NCBI rate limits:
      - Without API key: ~3 req/s  -> delay >= 0.34 s
      - With API key:    ~10 req/s -> delay >= 0.10 s

    Returns a mapping of {accession: faa_path} for successful downloads.
    """
    os.makedirs(fasta_dir, exist_ok=True)
    if ncbi_api_key and delay > 0.10:
        delay = 0.10

    downloaded: dict[str, str] = {}
    for i, acc in enumerate(accessions):
        out = os.path.join(fasta_dir, f"{acc}.faa")
        if os.path.exists(out):
            print(f"  [{i+1}/{len(accessions)}] {acc}: already downloaded, skipping")
            downloaded[acc] = out
            continue

        print(f"  [{i+1}/{len(accessions)}] {acc}: downloading protein FASTA...")
        path = download_species_faa(acc, out, ncbi_api_key)
        if path:
            downloaded[acc] = path
            record_count = sum(1 for _ in SeqIO.parse(path, "fasta"))
            print(f"  [{i+1}/{len(accessions)}] {acc}: saved {record_count} proteins -> {out}")

        if i < len(accessions) - 1:
            time.sleep(delay)

    print(f"\nDownloaded {len(downloaded)}/{len(accessions)} assemblies")
    return downloaded


# ---------------------------------------------------------------------------
# Parse FASTAs
# ---------------------------------------------------------------------------

def parse_all_fastas(fasta_paths: dict[str, str]) -> list[tuple[str, str]]:
    """
    Parse all FASTA files and return a flat list of (seq_id, aa_sequence) tuples.
    A single FASTA file may contain multiple records.
    """
    sequences: list[tuple[str, str]] = []
    for acc, path in fasta_paths.items():
        for record in SeqIO.parse(path, "fasta"):
            sequences.append((str(record.id), str(record.seq)))
    print(f"Parsed {len(sequences)} total sequences from {len(fasta_paths)} files")
    return sequences


# ---------------------------------------------------------------------------
# Combine FASTAs
# ---------------------------------------------------------------------------

def save_combined_fasta(sequences: list[tuple[str, str]], output_file: str) -> str:
    """Write all amino acid sequences to a single combined FASTA file."""
    with open(output_file, "w") as f:
        for seq_id, aa_seq in sequences:
            f.write(f">{seq_id}\n{aa_seq}\n")
    print(f"Combined {len(sequences)} sequences into {output_file}")
    return output_file


# ---------------------------------------------------------------------------
# Foldseek helpers
# ---------------------------------------------------------------------------

def resolve_foldseek_bin(foldseek_path: str = None) -> str:
    """Resolve the foldseek binary path."""
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


def _run(cmd: list[str], description: str):
    """Run a subprocess command with logging."""
    print(f"\n{description}: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        if result.stdout:
            print(result.stdout[:1000])
        if result.stderr:
            print("Info:", result.stderr[:500])
    except subprocess.CalledProcessError as e:
        print(f"Command failed:\n  stdout: {e.stdout}\n  stderr: {e.stderr}")
        raise


def download_prostt5_weights(
    foldseek_bin: str,
    weights_dir: str,
) -> str:
    """
    Download ProstT5 weights using foldseek's built-in downloader.

    Returns the path to the weights directory.
    """
    os.makedirs(weights_dir, exist_ok=True)
    weights_prefix = os.path.join(weights_dir, "prostt5")
    tmp_dir = os.path.join(weights_dir, "tmp_download")

    if os.path.exists(weights_prefix):
        print(f"ProstT5 weights already present at {weights_prefix}")
        return weights_prefix

    _run(
        [foldseek_bin, "databases", "ProstT5", weights_prefix, tmp_dir],
        "Downloading ProstT5 weights",
    )
    return weights_prefix


# ---------------------------------------------------------------------------
# 3Di generation via foldseek
# ---------------------------------------------------------------------------

def generate_3di_with_foldseek(
    combined_fasta: str,
    db_dir: str,
    foldseek_bin: str = "foldseek",
    prostt5_weights: str = None,
    use_gpu: bool = False,
    threads: int = 1,
    prostt5_split_length: int = None,
) -> tuple[str, str]:
    """
    Use foldseek createdb with --prostt5-model to predict 3Di codes from
    amino acid sequences.

    Note: Databases created this way contain only predicted 3Di structural
    sequences without Cα coordinates. This supports search and clustering
    but does not enable TM-score or LDDT output.

    Args:
        combined_fasta: Path to combined amino acid FASTA.
        db_dir: Directory for foldseek database files.
        foldseek_bin: Path to foldseek binary.
        prostt5_weights: Path prefix for ProstT5 weights (from `foldseek databases ProstT5`).
        use_gpu: Accelerate inference with GPU (--gpu 1).
        threads: Number of threads.

    Returns:
        tuple: (db_prefix, three_di_fasta_path)
    """
    os.makedirs(db_dir, exist_ok=True)
    db_prefix = os.path.join(db_dir, "sequenceDB")
    three_di_fasta = os.path.join(db_dir, "all_sequences_3di.fasta")

    # Create foldseek database with ProstT5 3Di prediction
    cmd = [
        foldseek_bin, "createdb",
        combined_fasta, db_prefix,
        "--threads", str(threads),
    ]
    if prostt5_weights:
        cmd.extend(["--prostt5-model", prostt5_weights])
    if use_gpu:
        cmd.extend(["--gpu", "1"])
    if prostt5_split_length is not None:
        cmd.extend(["-prostt5SplitLength", str(prostt5_split_length)])

    _run(cmd, "Creating foldseek database with 3Di prediction")

    # Link the sequence headers to the _ss (3Di) database so convert2fasta
    # produces FASTA records with the original sequence IDs
    header_db = db_prefix + "_h"
    ss_db = db_prefix + "_ss"
    ss_header = ss_db + "_h"

    _run(
        [foldseek_bin, "lndb", header_db, ss_header],
        "Linking headers to 3Di database",
    )

    # Extract 3Di sequences to FASTA format
    _run(
        [foldseek_bin, "convert2fasta", ss_db, three_di_fasta],
        "Extracting 3Di codes to FASTA",
    )

    # Count extracted sequences
    n = sum(1 for line in open(three_di_fasta) if line.startswith(">"))
    print(f"Saved {n} 3Di sequences to {three_di_fasta}")

    return db_prefix, three_di_fasta


# ---------------------------------------------------------------------------
# Foldseek all-vs-all
# ---------------------------------------------------------------------------

def run_foldseek_all_vs_all(
    db_prefix: str,
    output_dir: str,
    foldseek_bin: str = "foldseek",
    threads: int = 1,
    extra_args: list[str] = None,
) -> str:
    """
    Run foldseek search with the database against itself (all-vs-all),
    then convert results to a TSV file.

    Note: --alignment-type 1 (TMalign) is NOT supported for ProstT5-predicted
    databases (no Cα coordinates). The default alignment (3Di Gotoh-Smith-
    Waterman) is used instead.

    Returns the path to the result TSV file.
    """
    os.makedirs(output_dir, exist_ok=True)
    result_prefix = os.path.join(output_dir, "result")
    result_file = os.path.join(output_dir, "all_vs_all_results.tsv")
    tmp_dir = os.path.join(output_dir, "tmp_search")

    # All-vs-all search
    cmd = [
        foldseek_bin, "search",
        db_prefix, db_prefix,
        result_prefix, tmp_dir,
        "--threads", str(threads),
    ]
    if extra_args:
        cmd.extend(extra_args)

    _run(cmd, "Running foldseek all-vs-all search")

    # Convert to human-readable TSV
    _run(
        [foldseek_bin, "convertalis",
         db_prefix, db_prefix,
         result_prefix, result_file],
        "Converting results to TSV",
    )

    print(f"Results saved to {result_file}")
    return result_file


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Batch download whole-species protein FASTAs (.faa) from NCBI Datasets, "
            "convert to 3Di codes with foldseek (integrated ProstT5), then run "
            "all-vs-all foldseek search."
        )
    )
    parser.add_argument(
        "accession_file",
        help="Text file with one NCBI genome assembly accession (GCF_* or GCA_*) per line",
    )
    parser.add_argument(
        "--output-dir", "-o",
        default="batch_output",
        help="Directory for all output files (default: batch_output)",
    )
    parser.add_argument(
        "--fasta-dir",
        default=None,
        help="Directory containing pre-downloaded .faa files (skips NCBI download step)",
    )
    parser.add_argument(
        "--ncbi-api-key",
        default=None,
        help="NCBI API key for higher rate limits (optional)",
    )
    parser.add_argument(
        "--foldseek-path",
        default=None,
        help="Path to foldseek binary or installation directory",
    )
    parser.add_argument(
        "--prostt5-weights",
        default=None,
        help=(
            "Path prefix for pre-downloaded ProstT5 weights. "
            "If not provided, weights are downloaded automatically via "
            "'foldseek databases ProstT5' into the output directory."
        ),
    )
    parser.add_argument(
        "--gpu",
        action="store_true",
        help="Use GPU for ProstT5 inference (--gpu 1)",
    )
    parser.add_argument(
        "--prostt5-split-length",
        type=int, default=None,
        help="Max sequence length before splitting into multiple runs "
             "(-prostt5SplitLength). Shorter values reduce peak VRAM usage.",
    )
    parser.add_argument(
        "--threads", "-t",
        type=int,
        default=1,
        help="Number of threads for foldseek (default: 1)",
    )
    parser.add_argument(
        "--skip-foldseek",
        action="store_true",
        help="Stop after generating 3Di codes, skip foldseek search",
    )
    parser.add_argument(
        "--foldseek-args",
        nargs=argparse.REMAINDER,
        default=[],
        help="Extra arguments forwarded to foldseek search (e.g. --foldseek-args -e 0.001)",
    )
    args = parser.parse_args()

    out = Path(args.output_dir)
    fasta_dir = out / "fastas"
    combined_fasta = out / "combined_aa.fasta"
    db_dir = out / "foldseek_db"
    metadata_file = out / "metadata.json"
    out.mkdir(parents=True, exist_ok=True)

    # Resolve foldseek binary
    foldseek_bin = resolve_foldseek_bin(args.foldseek_path)

    # 1. Read accessions
    accessions = read_accessions(args.accession_file)

    # 2. Download FASTAs (or use pre-downloaded ones)
    if args.fasta_dir:
        print(f"\n--- Using pre-downloaded FASTAs from {args.fasta_dir} ---")
        fasta_dir = Path(args.fasta_dir)
        downloaded = {}
        for acc in accessions:
            faa_path = fasta_dir / f"{acc}.faa"
            if faa_path.exists():
                downloaded[acc] = str(faa_path)
            else:
                print(f"  [WARN] {acc}: no .faa file found at {faa_path}")
        print(f"Found {len(downloaded)}/{len(accessions)} pre-downloaded FASTAs")
    else:
        print("\n--- Downloading FASTAs ---")
        downloaded = batch_download(accessions, str(fasta_dir), args.ncbi_api_key)
    if not downloaded:
        print("No sequences found. Exiting.")
        return

    # 3. Parse
    print("\n--- Parsing sequences ---")
    sequences = parse_all_fastas(downloaded)
    if not sequences:
        print("No sequences parsed. Exiting.")
        return

    # 4. Combine into a single FASTA
    print("\n--- Combining sequences ---")
    save_combined_fasta(sequences, str(combined_fasta))

    # 5. Download ProstT5 weights if not provided
    prostt5_weights = args.prostt5_weights
    if not prostt5_weights:
        print("\n--- Downloading ProstT5 weights ---")
        prostt5_weights = download_prostt5_weights(
            foldseek_bin, str(out / "prostt5_weights"),
        )

    # 6. Generate 3Di codes via foldseek createdb
    print("\n--- Generating 3Di codes via foldseek ---")
    db_prefix, three_di_fasta = generate_3di_with_foldseek(
        str(combined_fasta),
        str(db_dir),
        foldseek_bin=foldseek_bin,
        prostt5_weights=prostt5_weights,
        use_gpu=args.gpu,
        threads=args.threads,
        prostt5_split_length=args.prostt5_split_length,
    )

    # 7. Foldseek all-vs-all
    result_file = None
    if args.skip_foldseek:
        print("\nSkipping foldseek search (--skip-foldseek set).")
    else:
        print("\n--- Running foldseek all-vs-all ---")
        try:
            result_file = run_foldseek_all_vs_all(
                db_prefix,
                output_dir=str(out / "foldseek_results"),
                foldseek_bin=foldseek_bin,
                threads=args.threads,
                extra_args=args.foldseek_args or None,
            )
        except (FileNotFoundError, subprocess.CalledProcessError):
            print(f"\nFoldseek search failed. 3Di sequences are available at {three_di_fasta}")

    # 8. Save metadata
    metadata = {
        "accession_file": args.accession_file,
        "assemblies_requested": accessions,
        "assemblies_downloaded": list(downloaded.keys()),
        "num_sequences": len(sequences),
        "three_di_fasta": str(three_di_fasta),
        "foldseek_results": result_file,
    }
    with open(metadata_file, "w") as f:
        json.dump(metadata, f, indent=2)
    print(f"\nDone. All outputs in {out}/")


if __name__ == "__main__":
    main()
