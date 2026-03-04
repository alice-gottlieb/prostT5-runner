"""
Batch download whole-species protein FASTA (.faa) files from NCBI Datasets,
convert all sequences to 3Di codes with ProstT5, then run an all-vs-all
foldseek search across all sequences.

Input:  a text file with one NCBI genome assembly accession (GCF_* or GCA_*)
        per line
Output: per-assembly .faa files, a combined 3Di FASTA, and foldseek results
"""

import re
import io
import time
import json
import os
import subprocess
import argparse
import zipfile
import requests
import torch
from pathlib import Path
from Bio import SeqIO
from transformers import T5Tokenizer, AutoModelForSeq2SeqLM


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
# 3Di generation
# ---------------------------------------------------------------------------

def generate_3di_codes(
    sequences: list[tuple[str, str]],
    device: str = None,
    batch_size: int = 1,
    model_dir: str = None,
) -> dict[str, str]:
    """
    Convert amino acid sequences to 3Di codes using ProstT5.

    Loads the model once and processes sequences one at a time (ProstT5
    requires per-sequence length constraints for min/max_length).

    Args:
        model_dir: Local directory to load weights from (if it already contains
                   a saved model) or save weights to after downloading. If None,
                   weights are downloaded fresh each run into HuggingFace's
                   default cache (~/.cache/huggingface/).

    Returns a dict mapping seq_id -> 3Di string.
    """
    if device is None:
        device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    else:
        device = torch.device(device)

    print(f"Using device: {device}")

    # Determine whether to load from a local directory or download from HuggingFace
    model_source = "Rostlab/ProstT5"
    if model_dir and os.path.isfile(os.path.join(model_dir, "config.json")):
        model_source = model_dir
        print(f"Loading ProstT5 from local directory: {model_dir}")
    else:
        print("Loading ProstT5 model and tokenizer from HuggingFace...")

    tokenizer = T5Tokenizer.from_pretrained(model_source, do_lower_case=False)
    model = AutoModelForSeq2SeqLM.from_pretrained(model_source).to(device)

    # Save weights locally for future runs if a directory was specified but didn't exist yet
    if model_dir and model_source != model_dir:
        os.makedirs(model_dir, exist_ok=True)
        print(f"Saving model weights to {model_dir} for future runs...")
        tokenizer.save_pretrained(model_dir)
        model.save_pretrained(model_dir)
    model.float()
    model.eval()

    gen_kwargs = {
        "do_sample": True,
        "num_beams": 3,
        "top_p": 0.95,
        "temperature": 1.2,
        "top_k": 6,
        "repetition_penalty": 1.2,
    }

    results: dict[str, str] = {}
    print(f"Generating 3Di codes for {len(sequences)} sequences...")

    for i, (seq_id, aa_seq) in enumerate(sequences):
        aa_seq = re.sub(r"[UZOB]", "X", aa_seq.upper())
        input_seq = "<AA2fold> " + " ".join(list(aa_seq))
        seq_len = len(aa_seq)

        ids = tokenizer(
            [input_seq],
            add_special_tokens=True,
            padding="longest",
            return_tensors="pt",
        ).to(device)

        with torch.no_grad():
            out = model.generate(
                ids.input_ids,
                attention_mask=ids.attention_mask,
                max_length=seq_len,
                min_length=seq_len,
                **gen_kwargs,
            )

        three_di = tokenizer.batch_decode(out, skip_special_tokens=True)[0]
        three_di = three_di.replace(" ", "").lower()
        results[seq_id] = three_di
        print(f"  [{i+1}/{len(sequences)}] {seq_id}: {seq_len} AA -> {len(three_di)} 3Di")

    return results


# ---------------------------------------------------------------------------
# Save 3Di
# ---------------------------------------------------------------------------

def save_3di_fasta(three_di_codes: dict[str, str], output_file: str) -> str:
    """Write all 3Di codes to a single combined FASTA file."""
    with open(output_file, "w") as f:
        for seq_id, code in three_di_codes.items():
            f.write(f">{seq_id}\n{code}\n")
    print(f"Saved {len(three_di_codes)} 3Di sequences to {output_file}")
    return output_file


# ---------------------------------------------------------------------------
# Foldseek all-vs-all
# ---------------------------------------------------------------------------

def run_foldseek_all_vs_all(
    three_di_fasta: str,
    output_dir: str,
    foldseek_path: str = None,
    tmp_dir: str = None,
    extra_args: list[str] = None,
) -> str:
    """
    Run foldseek easy-search with the 3Di FASTA as both query and target,
    producing an all-vs-all structural similarity comparison.

    Returns the path to the result file.
    """
    os.makedirs(output_dir, exist_ok=True)
    result_file = os.path.join(output_dir, "all_vs_all_results.tsv")
    tmp = tmp_dir or os.path.join(output_dir, "foldseek_tmp")

    # Resolve foldseek binary
    if foldseek_path:
        if os.path.isdir(foldseek_path):
            bin_path = os.path.join(foldseek_path, "bin", "foldseek")
            if not os.path.exists(bin_path):
                bin_path = os.path.join(foldseek_path, "foldseek")
        else:
            bin_path = foldseek_path
        if not os.path.exists(bin_path):
            raise FileNotFoundError(f"foldseek binary not found at {bin_path}")
    else:
        bin_path = "foldseek"

    cmd = [
        bin_path, "easy-search",
        three_di_fasta,   # query
        three_di_fasta,   # target (same file = all-vs-all)
        result_file,
        tmp,
        "--alignment-type", "1",   # 1 = 3Di alignment
    ]
    if extra_args:
        cmd.extend(extra_args)

    print(f"\nRunning foldseek: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print("Foldseek completed successfully.")
        if result.stdout:
            print(result.stdout)
        if result.stderr:
            print("Info:", result.stderr[:500])
    except subprocess.CalledProcessError as e:
        print(f"Foldseek failed:\n{e.stderr}")
        raise
    except FileNotFoundError:
        print("foldseek not found. Install it or provide --foldseek-path.")
        raise

    print(f"Results saved to {result_file}")
    return result_file


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Batch download whole-species protein FASTAs (.faa) from NCBI Datasets, "
            "convert to 3Di codes with ProstT5, then run all-vs-all foldseek search."
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
        "--ncbi-api-key",
        default=None,
        help="NCBI API key for higher rate limits (optional)",
    )
    parser.add_argument(
        "--device",
        choices=["cuda", "cpu"],
        default=None,
        help="Device for ProstT5 inference (auto-detected if not set)",
    )
    parser.add_argument(
        "--model-dir",
        default=None,
        help=(
            "Directory to save/load ProstT5 weights. "
            "First run downloads and saves here; subsequent runs load from disk."
        ),
    )
    parser.add_argument(
        "--foldseek-path",
        default=None,
        help="Path to foldseek binary or installation directory",
    )
    parser.add_argument(
        "--skip-foldseek",
        action="store_true",
        help="Stop after generating 3Di codes, skip foldseek",
    )
    parser.add_argument(
        "--foldseek-args",
        nargs=argparse.REMAINDER,
        default=[],
        help="Extra arguments forwarded to foldseek (e.g. --foldseek-args -e 0.001)",
    )
    args = parser.parse_args()

    out = Path(args.output_dir)
    fasta_dir = out / "fastas"
    three_di_fasta = out / "all_sequences_3di.fasta"
    metadata_file = out / "metadata.json"
    out.mkdir(parents=True, exist_ok=True)

    # 1. Read accessions
    accessions = read_accessions(args.accession_file)

    # 2. Download FASTAs
    print("\n--- Downloading FASTAs ---")
    downloaded = batch_download(accessions, str(fasta_dir), args.ncbi_api_key)
    if not downloaded:
        print("No sequences downloaded. Exiting.")
        return

    # 3. Parse
    print("\n--- Parsing sequences ---")
    sequences = parse_all_fastas(downloaded)
    if not sequences:
        print("No sequences parsed. Exiting.")
        return

    # 4. Generate 3Di codes
    print("\n--- Generating 3Di codes ---")
    three_di_codes = generate_3di_codes(sequences, device=args.device, model_dir=args.model_dir)

    # 5. Save combined 3Di FASTA
    print("\n--- Saving 3Di FASTA ---")
    save_3di_fasta(three_di_codes, str(three_di_fasta))

    # 6. Foldseek all-vs-all
    result_file = None
    if args.skip_foldseek:
        print("\nSkipping foldseek (--skip-foldseek set).")
    else:
        print("\n--- Running foldseek all-vs-all ---")
        try:
            result_file = run_foldseek_all_vs_all(
                str(three_di_fasta),
                output_dir=str(out / "foldseek_results"),
                foldseek_path=args.foldseek_path,
                extra_args=args.foldseek_args or None,
            )
        except (FileNotFoundError, subprocess.CalledProcessError):
            print(f"\nFoldseek failed. 3Di sequences are available at {three_di_fasta}")
            print("You can run manually:")
            print(f"  foldseek easy-search {three_di_fasta} {three_di_fasta} results.tsv tmp --alignment-type 1")

    # 7. Save metadata
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
