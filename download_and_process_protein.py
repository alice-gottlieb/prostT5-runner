"""
Download a bacterial protein FASTA from NCBI datasets and process with ProstT5 and foldseek.

This script:
1. Downloads a bacterial protein FASTA file from NCBI datasets
2. Generates 3Di code using ProstT5 model
3. Runs the 3Di sequences through foldseek for structure search
"""

import requests
import re
import torch
import subprocess
import tempfile
import os
from pathlib import Path
from Bio import SeqIO
from transformers import T5Tokenizer, AutoModelForSeq2SeqLM
import argparse
import json


def download_fasta_from_ncbi(accession_id: str, output_file: str = None) -> str:
    """
    Download a protein FASTA sequence from NCBI datasets.
    
    Args:
        accession_id: The accession ID (e.g., 'NP_000001' or 'GCF_000001405.40')
        output_file: Optional output file path. If None, creates temp file.
    
    Returns:
        Path to the downloaded FASTA file
    """
    print(f"Downloading protein sequence for {accession_id} from NCBI...")
    
    # Use NCBI Entrez API to download the sequence
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "protein",
        "id": accession_id,
        "rettype": "fasta",
        "retmode": "text"
    }
    
    response = requests.get(base_url, params=params)
    response.raise_for_status()
    
    if output_file is None:
        output_file = f"protein_{accession_id}.fasta"
    
    with open(output_file, 'w') as f:
        f.write(response.text)
    
    print(f"Downloaded FASTA to {output_file}")
    return output_file


def parse_fasta(fasta_file: str) -> list:
    """
    Parse a FASTA file and extract sequences.
    
    Args:
        fasta_file: Path to the FASTA file
    
    Returns:
        List of tuples (sequence_id, sequence)
    """
    sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append((str(record.id), str(record.seq)))
    
    print(f"Parsed {len(sequences)} sequences from {fasta_file}")
    return sequences


def generate_3di_codes(sequences: list, device: str = None) -> dict:
    """
    Generate 3Di codes from amino acid sequences using ProstT5.
    
    Args:
        sequences: List of tuples (seq_id, aa_sequence)
        device: Device to use ('cuda' or 'cpu'). Auto-detects if None.
    
    Returns:
        Dictionary mapping sequence IDs to 3Di codes
    """
    if device is None:
        device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    else:
        device = torch.device(device)
    
    print(f"Using device: {device}")
    print("Loading ProstT5 model and tokenizer...")
    
    # Load the tokenizer and model
    tokenizer = T5Tokenizer.from_pretrained('Rostlab/ProstT5', do_lower_case=False)
    model = AutoModelForSeq2SeqLM.from_pretrained("Rostlab/ProstT5").to(device)
    
    # Use mixed precision if GPU available
    model.float() if device.type == 'cpu' else model.half()
    
    print("Generating 3Di codes...")
    
    results = {}
    
    # Process in batches
    for seq_id, aa_sequence in sequences:
        # Replace rare/ambiguous amino acids
        aa_sequence = re.sub(r"[UZOB]", "X", aa_sequence.upper())
        
        # Add spaces between amino acids
        spaced_sequence = " ".join(list(aa_sequence))
        
        # Add prefix for AA to 3Di translation
        input_sequence = "<AA2fold> " + spaced_sequence
        
        min_len = len(aa_sequence)
        max_len = len(aa_sequence)
        
        # Tokenize
        ids = tokenizer([input_sequence],
                       add_special_tokens=True,
                       padding="longest",
                       return_tensors='pt').to(device)
        
        # Generation configuration
        gen_kwargs = {
            "do_sample": True,
            "num_beams": 3,
            "top_p": 0.95,
            "temperature": 1.2,
            "top_k": 6,
            "repetition_penalty": 1.2,
        }
        
        # Generate 3Di code
        with torch.no_grad():
            translations = model.generate(
                ids.input_ids,
                attention_mask=ids.attention_mask,
                max_length=max_len,
                min_length=min_len,
                **gen_kwargs
            )
        
        # Decode the translation
        three_di_code = tokenizer.batch_decode(translations, skip_special_tokens=True)
        three_di_code = three_di_code[0].replace(" ", "").lower()
        
        results[seq_id] = three_di_code
        print(f"  {seq_id}: {len(aa_sequence)} AAs -> {len(three_di_code)} 3Di chars")
    
    return results


def save_3di_fasta(three_di_codes: dict, output_file: str = "sequences_3di.fasta") -> str:
    """
    Save 3Di codes to a FASTA file.
    
    Args:
        three_di_codes: Dictionary mapping sequence IDs to 3Di codes
        output_file: Output FASTA file path
    
    Returns:
        Path to the output file
    """
    with open(output_file, 'w') as f:
        for seq_id, three_di_code in three_di_codes.items():
            f.write(f">{seq_id}\n{three_di_code}\n")
    
    print(f"Saved 3Di sequences to {output_file}")
    return output_file


def run_foldseek(query_file: str, target_db: str = None, output_dir: str = None) -> str:
    """
    Run foldseek to search 3Di sequences against a database.
    
    Args:
        query_file: Path to the query 3Di FASTA file
        target_db: Path to the target database. If None, searches against environment default.
        output_dir: Output directory for foldseek results
    
    Returns:
        Path to the results
    """
    if output_dir is None:
        output_dir = "foldseek_results"
    
    os.makedirs(output_dir, exist_ok=True)
    result_file = os.path.join(output_dir, "search_results")
    
    # Build foldseek command
    cmd = ["foldseek", "easy-search", query_file]
    
    if target_db:
        cmd.append(target_db)
    else:
        # Search against PDB database (requires foldseek setup)
        cmd.append("PDB")
    
    cmd.extend([result_file, "/tmp/foldseek_tmp"])
    
    print(f"Running foldseek command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print("Foldseek completed successfully!")
        print(result.stdout)
        if result.stderr:
            print("Warnings/Info:", result.stderr)
    except subprocess.CalledProcessError as e:
        print(f"Error running foldseek: {e.stderr}")
        raise
    except FileNotFoundError:
        print("Error: foldseek is not installed or not in PATH.")
        print("Please install foldseek: https://github.com/steineggerlab/foldseek")
        raise
    
    # List output files
    output_files = [f for f in os.listdir(output_dir) if "search_results" in f]
    if output_files:
        print(f"Results saved to: {output_dir}/")
        for f in output_files:
            print(f"  - {f}")
    
    return result_file


def main():
    parser = argparse.ArgumentParser(
        description="Download bacterial protein FASTA and process with ProstT5 and foldseek"
    )
    parser.add_argument(
        "--accession",
        type=str,
        help="NCBI protein accession ID (e.g., 'NP_000001')"
    )
    parser.add_argument(
        "--fasta",
        type=str,
        help="Path to local FASTA file (if not downloading)"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="protein_processing_output",
        help="Output directory for results"
    )
    parser.add_argument(
        "--foldseek-db",
        type=str,
        default=None,
        help="Foldseek database path or name"
    )
    parser.add_argument(
        "--device",
        type=str,
        choices=['cuda', 'cpu'],
        default=None,
        help="Device to use for ProstT5 (auto-detected if not specified)"
    )
    parser.add_argument(
        "--skip-foldseek",
        action="store_true",
        help="Skip foldseek step"
    )
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Step 1: Get FASTA file
    if args.accession:
        fasta_file = download_fasta_from_ncbi(
            args.accession,
            os.path.join(args.output_dir, f"protein_{args.accession}.fasta")
        )
    elif args.fasta:
        fasta_file = args.fasta
        print(f"Using local FASTA file: {fasta_file}")
    else:
        parser.error("Either --accession or --fasta must be provided")
    
    # Step 2: Parse FASTA
    sequences = parse_fasta(fasta_file)
    
    # Step 3: Generate 3Di codes
    three_di_codes = generate_3di_codes(sequences, device=args.device)
    
    # Step 4: Save 3Di FASTA
    three_di_fasta = save_3di_fasta(
        three_di_codes,
        os.path.join(args.output_dir, "sequences_3di.fasta")
    )
    
    # Step 5: Run foldseek (optional)
    if not args.skip_foldseek:
        run_foldseek(three_di_fasta, target_db=args.foldseek_db, output_dir=args.output_dir)
    
    # Save metadata
    metadata = {
        "input_fasta": fasta_file,
        "num_sequences": len(sequences),
        "three_di_fasta": three_di_fasta,
    }
    
    with open(os.path.join(args.output_dir, "metadata.json"), 'w') as f:
        json.dump(metadata, f, indent=2)
    
    print(f"\nProcessing complete! Results saved to {args.output_dir}/")


if __name__ == "__main__":
    main()
