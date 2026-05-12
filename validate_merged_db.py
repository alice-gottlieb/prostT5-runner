"""
Validate that a merged foldseek DB contains every protein from a folder of
.faa files and that each 3Di code has the same length as its source protein.

Checks performed for every protein ID found in the .faa files:
  1. MISSING   — the ID is not present in the merged DB's 3Di sequences.
  2. LENGTH    — the 3Di sequence length differs from the amino-acid length.

Output is written to a log file (default: <db_dir>/validation.log) and
a summary is printed to the terminal.

Usage
-----
    python validate_merged_db.py --db-dir /scratch/mergedDB_dir \\
                                 --faa-dir /scratch/genomes \\
                                 [--foldseek  /path/to/foldseek] \\
                                 [--log-file  /path/to/validation.log]
"""

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

from Bio import SeqIO


def load_faa_dir(faa_dir: Path) -> dict[str, int]:
    """Read all .faa files in faa_dir and return {protein_id: aa_length}."""
    proteins: dict[str, int] = {}
    faa_files = sorted(faa_dir.glob("*.faa"))
    if not faa_files:
        sys.exit(f"No .faa files found in {faa_dir}")
    for faa_path in faa_files:
        for record in SeqIO.parse(faa_path, "fasta"):
            proteins[record.id] = len(record.seq)
    return proteins


def extract_3di(db_dir: Path, foldseek_bin: str) -> dict[str, int]:
    """Load the merged 3Di FASTA, generating it from the sub-DB if needed.
    {protein_id: 3di_length}."""
    precomputed_fasta = db_dir / "all_sequences_3di.fasta"
    if precomputed_fasta.exists():
        print(f"  Using existing 3Di FASTA: {precomputed_fasta}")
        return {
            record.id: len(record.seq)
            for record in SeqIO.parse(precomputed_fasta, "fasta")
        }

    ss_db = db_dir / "mergedDB_ss"
    if not ss_db.exists():
        sys.exit(f"3Di sub-DB not found: {ss_db}")

    with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as tmp:
        tmp_path = Path(tmp.name)

    try:
        subprocess.run(
            [foldseek_bin, "convert2fasta", str(ss_db), str(tmp_path)],
            check=True,
            capture_output=True,
            text=True,
        )
        di3: dict[str, int] = {}
        for record in SeqIO.parse(tmp_path, "fasta"):
            di3[record.id] = len(record.seq)
    finally:
        tmp_path.unlink(missing_ok=True)

    return di3


def validate(proteins: dict[str, int], di3: dict[str, int]) -> tuple[list, list]:
    """Return (missing_ids, length_mismatches) where:
      missing_ids:       [protein_id, ...]
      length_mismatches: [(protein_id, aa_len, di3_len), ...]
    """
    missing = []
    mismatches = []
    for prot_id, aa_len in proteins.items():
        if prot_id not in di3:
            missing.append(prot_id)
        elif di3[prot_id] != aa_len:
            mismatches.append((prot_id, aa_len, di3[prot_id]))
    return missing, mismatches


def write_log(log_path: Path,
              proteins: dict[str, int],
              di3: dict[str, int],
              missing: list,
              mismatches: list) -> None:
    extra_in_db = sorted(set(di3) - set(proteins))
    with open(log_path, "w") as f:
        f.write("=== Validation summary ===\n")
        f.write(f"  Proteins in .faa files : {len(proteins)}\n")
        f.write(f"  Sequences in merged DB : {len(di3)}\n")
        f.write(f"  Missing from DB        : {len(missing)}\n")
        f.write(f"  Length mismatches      : {len(mismatches)}\n")
        f.write(f"  Extra in DB (not in faa): {len(extra_in_db)}\n")
        f.write("\n")

        if missing:
            f.write("=== MISSING proteins (in .faa but not in DB) ===\n")
            for prot_id in sorted(missing):
                f.write(f"  MISSING  {prot_id}  (aa_len={proteins[prot_id]})\n")
            f.write("\n")

        if mismatches:
            f.write("=== LENGTH MISMATCH proteins ===\n")
            for prot_id, aa_len, di3_len in sorted(mismatches):
                f.write(f"  MISMATCH {prot_id}  "
                        f"aa_len={aa_len}  3di_len={di3_len}  "
                        f"diff={di3_len - aa_len:+d}\n")
            f.write("\n")

        if extra_in_db:
            f.write("=== Extra sequences in DB (not in any .faa file) ===\n")
            for prot_id in extra_in_db:
                f.write(f"  EXTRA    {prot_id}\n")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--db-dir", required=True, type=Path,
                    help="Directory containing the merged foldseek DB "
                         "(must have mergedDB_ss and mergedDB_ss_h)")
    ap.add_argument("--faa-dir", required=True, type=Path,
                    help="Directory containing .faa protein files")
    ap.add_argument("--foldseek", default="foldseek",
                    help="Path to foldseek binary (default: foldseek on PATH)")
    ap.add_argument("--log-file", default=None, type=Path,
                    help="Path for the validation log "
                         "(default: <db_dir>/validation.log)")
    args = ap.parse_args()

    log_path = args.log_file or (args.db_dir / "validation.log")

    print(f"Loading proteins from {args.faa_dir} ...")
    proteins = load_faa_dir(args.faa_dir)
    print(f"  {len(proteins)} proteins loaded")

    print(f"Extracting 3Di sequences from {args.db_dir} ...")
    di3 = extract_3di(args.db_dir, args.foldseek)
    print(f"  {len(di3)} 3Di sequences extracted")

    print("Validating ...")
    missing, mismatches = validate(proteins, di3)

    write_log(log_path, proteins, di3, missing, mismatches)

    extra_in_db = len(set(di3) - set(proteins))
    print()
    print("=== Results ===")
    print(f"  Proteins in .faa files : {len(proteins)}")
    print(f"  Sequences in merged DB : {len(di3)}")
    print(f"  Missing from DB        : {len(missing)}")
    print(f"  Length mismatches      : {len(mismatches)}")
    print(f"  Extra in DB (not in faa): {extra_in_db}")
    print()
    print(f"Full details written to: {log_path}")

    # if missing or mismatches:
    #     sys.exit(1)


if __name__ == "__main__":
    main()
