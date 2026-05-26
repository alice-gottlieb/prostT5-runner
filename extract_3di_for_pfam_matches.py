#!/usr/bin/env python3
"""Add matching 3Di subsequences to HMMER hmmsearch domtblout files.

For each non-comment row in a --domtblout file, the script uses the target
sequence ID and HMMER target-sequence coordinates to extract the corresponding
region from a 3Di FASTA file. By default it extracts the aligned coordinates
("ali coord"); use --coords env to extract the domain envelope instead.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable


MIN_DOMTBLOUT_FIELDS = 22


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Append a 3di_sequence column to HMMER hmmsearch --domtblout files."
        )
    )
    parser.add_argument(
        "domtblout",
        type=Path,
        help="Input .domtblout file, or a directory containing *.domtblout files.",
    )
    parser.add_argument(
        "three_di_fasta",
        type=Path,
        help="FASTA file containing full-length 3Di sequences keyed by protein ID.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help=(
            "Output file for a single domtblout input. Defaults to "
            "<input>.with_3di.domtblout."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help=(
            "Output directory for a domtblout directory input. Defaults to the "
            "input directory."
        ),
    )
    parser.add_argument(
        "--coords",
        choices=("ali", "env"),
        default="ali",
        help="Target-sequence coordinates to extract: aligned domain or envelope.",
    )
    parser.add_argument(
        "--missing-value",
        default="",
        help="Value written when a target ID is absent from the 3Di FASTA.",
    )
    return parser.parse_args()


def find_domtblout_files(path: Path) -> list[Path]:
    if path.is_dir():
        return sorted(path.glob("*.domtblout"))
    return [path]


def output_path_for(input_path: Path, args: argparse.Namespace) -> Path:
    if args.output is not None:
        return args.output

    stem = input_path.name
    if stem.endswith(".domtblout"):
        stem = stem[: -len(".domtblout")]
    else:
        stem = input_path.stem

    output_dir = args.output_dir if args.output_dir is not None else input_path.parent
    return output_dir / f"{stem}.with_3di.domtblout"


def iter_data_fields(domtblout: Path) -> Iterable[list[str]]:
    with domtblout.open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split(maxsplit=MIN_DOMTBLOUT_FIELDS)
            if len(fields) < MIN_DOMTBLOUT_FIELDS:
                continue
            yield fields


def collect_needed_ids(domtblout_files: Iterable[Path]) -> set[str]:
    needed: set[str] = set()
    for domtblout in domtblout_files:
        for fields in iter_data_fields(domtblout):
            needed.add(fields[0])
    return needed


def read_needed_fasta_sequences(fasta_path: Path, needed_ids: set[str]) -> dict[str, str]:
    sequences: dict[str, str] = {}
    current_id: str | None = None
    current_parts: list[str] = []
    keep_current = False

    def flush_current() -> None:
        if current_id is not None and keep_current:
            sequences[current_id] = "".join(current_parts)

    with fasta_path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush_current()
                current_id = line[1:].split()[0]
                current_parts = []
                keep_current = current_id in needed_ids
            elif keep_current:
                current_parts.append(line)

    flush_current()
    return sequences


def extract_subsequence(sequence: str, start: int, end: int) -> str:
    if start > end:
        start, end = end, start
    start = max(start, 1)
    end = min(end, len(sequence))
    if start > end:
        return ""
    return sequence[start - 1 : end]


def coords_from_fields(fields: list[str], coord_set: str) -> tuple[int, int]:
    if coord_set == "ali":
        return int(fields[17]), int(fields[18])
    return int(fields[19]), int(fields[20])


def append_3di_column(
    domtblout: Path,
    output_path: Path,
    sequences: dict[str, str],
    coord_set: str,
    missing_value: str,
) -> tuple[int, int]:
    rows_written = 0
    missing_ids = 0

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with domtblout.open() as in_handle, output_path.open("w") as out_handle:
        for line in in_handle:
            stripped = line.rstrip("\n")
            if line.startswith("#"):
                if stripped.startswith("# target name"):
                    out_handle.write(f"{stripped}\t3di_sequence\n")
                else:
                    out_handle.write(line)
                continue

            if not stripped:
                out_handle.write(line)
                continue

            fields = stripped.split(maxsplit=MIN_DOMTBLOUT_FIELDS)
            if len(fields) < MIN_DOMTBLOUT_FIELDS:
                out_handle.write(f"{stripped}\t{missing_value}\n")
                continue

            target_id = fields[0]
            sequence = sequences.get(target_id)
            if sequence is None:
                subsequence = missing_value
                missing_ids += 1
            else:
                start, end = coords_from_fields(fields, coord_set)
                subsequence = extract_subsequence(sequence, start, end)

            out_handle.write(f"{stripped}\t{subsequence}\n")
            rows_written += 1

    return rows_written, missing_ids


def main() -> int:
    args = parse_args()

    domtblout_files = find_domtblout_files(args.domtblout)
    if not domtblout_files:
        raise SystemExit(f"No .domtblout files found: {args.domtblout}")
    if args.domtblout.is_dir() and args.output is not None:
        raise SystemExit("--output can only be used with a single domtblout file")

    for path in domtblout_files:
        if not path.is_file():
            raise SystemExit(f"domtblout file not found: {path}")
    if not args.three_di_fasta.is_file():
        raise SystemExit(f"3Di FASTA file not found: {args.three_di_fasta}")

    needed_ids = collect_needed_ids(domtblout_files)
    print(f"Found {len(needed_ids)} unique target IDs in {len(domtblout_files)} domtblout file(s).")
    sequences = read_needed_fasta_sequences(args.three_di_fasta, needed_ids)
    print(f"Loaded {len(sequences)} matching 3Di sequence(s) from {args.three_di_fasta}.")

    total_rows = 0
    total_missing = 0
    for domtblout in domtblout_files:
        output_path = output_path_for(domtblout, args)
        rows, missing = append_3di_column(
            domtblout=domtblout,
            output_path=output_path,
            sequences=sequences,
            coord_set=args.coords,
            missing_value=args.missing_value,
        )
        total_rows += rows
        total_missing += missing
        print(f"Wrote {rows} row(s) to {output_path} ({missing} missing 3Di IDs).")

    if total_missing:
        print(f"WARN: {total_missing} domain row(s) did not have a matching 3Di sequence.")
    print(f"Done. Added 3Di subsequences for {total_rows - total_missing}/{total_rows} domain row(s).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
