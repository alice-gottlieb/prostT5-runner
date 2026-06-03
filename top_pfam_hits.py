#!/usr/bin/env python3
"""Find the top HMMER hits for a Pfam domain in pfam_hits2-1by1.

Examples:
    uv run python top_pfam_hits.py AAA 10
    uv run python top_pfam_hits.py PF00004 10 --unique-target
    uv run python top_pfam_hits.py PF00004.36 10 --score full
    uv run python top_pfam_hits.py PF00198 5 -i pfam_hits2-1by1/GCF_000005845.2.with_3di.domtblout
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


MIN_DOMTBLOUT_FIELDS = 22


@dataclass(frozen=True)
class PfamHit:
    source_file: str
    target_name: str
    target_accession: str
    target_length: int
    domain_name: str
    domain_accession: str
    domain_length: int
    full_evalue: float
    full_score: float
    full_bias: float
    domain_number: int
    domain_count: int
    conditional_evalue: float
    independent_evalue: float
    domain_score: float
    domain_bias: float
    hmm_from: int
    hmm_to: int
    ali_from: int
    ali_to: int
    env_from: int
    env_to: int
    accuracy: float
    description: str
    three_di_sequence: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Find the top N rows for a Pfam domain across HMMER domtblout files "
            "in pfam_hits2-1by1."
        )
    )
    parser.add_argument(
        "domain",
        help=(
            "Pfam domain query name or accession. Examples: AAA, PF00004, "
            "PF00004.36."
        ),
    )
    parser.add_argument("n", type=int, help="Number of top hits to print.")
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        default=Path("pfam_hits2-1by1"),
        help=(
            "A .with_3di.domtblout file, or a directory containing "
            "*.with_3di.domtblout files."
        ),
    )
    parser.add_argument(
        "--score",
        choices=("full", "domain"),
        default="domain",
        help=(
            "Rank by full-sequence score or per-domain score. Ties are broken "
            "by the corresponding e-value. Defaults to domain."
        ),
    )
    parser.add_argument(
        "--unique-target",
        action="store_true",
        help="Keep only the best row per target protein before selecting top N.",
    )
    parser.add_argument(
        "--include-plain",
        action="store_true",
        help="Also include plain *.domtblout files when the input is a directory.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Write TSV output to this path instead of stdout.",
    )
    return parser.parse_args()


def normalize_domain(value: str) -> str:
    value = value.strip().lower()
    if value.startswith("pf") and "." in value:
        return value.split(".", 1)[0]
    return value


def domain_matches(hit: PfamHit, requested_domain: str) -> bool:
    requested = normalize_domain(requested_domain)
    domain_name = normalize_domain(hit.domain_name)
    domain_accession = normalize_domain(hit.domain_accession)
    return requested in {domain_name, domain_accession}


def domtblout_files(input_path: Path, include_plain: bool) -> list[Path]:
    if input_path.is_file():
        return [input_path]
    if not input_path.is_dir():
        raise SystemExit(f"Input path not found: {input_path}")

    paths = sorted(input_path.glob("*.with_3di.domtblout"))
    if include_plain:
        plain_paths = [
            path
            for path in sorted(input_path.glob("*.domtblout"))
            if not path.name.endswith(".with_3di.domtblout")
        ]
        paths.extend(plain_paths)
    return paths


def split_description_and_3di(rest: str) -> tuple[str, str]:
    if "\t" not in rest:
        return rest, ""
    description, three_di_sequence = rest.rsplit("\t", 1)
    return description, three_di_sequence


def parse_domtblout_row(line: str, source_file: Path) -> PfamHit | None:
    fields = line.rstrip("\n").split(maxsplit=MIN_DOMTBLOUT_FIELDS)
    if len(fields) < MIN_DOMTBLOUT_FIELDS + 1:
        return None
    description, three_di_sequence = split_description_and_3di(fields[22])

    return PfamHit(
        source_file=source_file.name,
        target_name=fields[0],
        target_accession=fields[1],
        target_length=int(fields[2]),
        domain_name=fields[3],
        domain_accession=fields[4],
        domain_length=int(fields[5]),
        full_evalue=float(fields[6]),
        full_score=float(fields[7]),
        full_bias=float(fields[8]),
        domain_number=int(fields[9]),
        domain_count=int(fields[10]),
        conditional_evalue=float(fields[11]),
        independent_evalue=float(fields[12]),
        domain_score=float(fields[13]),
        domain_bias=float(fields[14]),
        hmm_from=int(fields[15]),
        hmm_to=int(fields[16]),
        ali_from=int(fields[17]),
        ali_to=int(fields[18]),
        env_from=int(fields[19]),
        env_to=int(fields[20]),
        accuracy=float(fields[21]),
        description=description,
        three_di_sequence=three_di_sequence,
    )


def iter_hits(paths: Iterable[Path]) -> Iterable[PfamHit]:
    for path in paths:
        with path.open() as handle:
            for line_number, line in enumerate(handle, start=1):
                if line.startswith("#") or not line.strip():
                    continue
                try:
                    hit = parse_domtblout_row(line, path)
                except ValueError as exc:
                    print(
                        f"WARN: skipped malformed row {path}:{line_number}: {exc}",
                        file=sys.stderr,
                    )
                    continue
                if hit is not None:
                    yield hit


def hit_sort_key(hit: PfamHit, score: str) -> tuple[float, float, str, int]:
    if score == "domain":
        return (-hit.domain_score, hit.independent_evalue, hit.target_name, hit.domain_number)
    return (-hit.full_score, hit.full_evalue, hit.target_name, hit.domain_number)


def best_per_target(hits: Iterable[PfamHit], score: str) -> list[PfamHit]:
    best: dict[str, PfamHit] = {}
    for hit in hits:
        current = best.get(hit.target_name)
        if current is None or hit_sort_key(hit, score) < hit_sort_key(current, score):
            best[hit.target_name] = hit
    return list(best.values())


def select_top_hits(hits: Iterable[PfamHit], score: str, n: int) -> list[PfamHit]:
    return sorted(hits, key=lambda hit: hit_sort_key(hit, score))[:n]


def write_hits(hits: list[PfamHit], output_path: Path | None) -> None:
    columns = [
        "rank",
        "target_name",
        "domain_name",
        "domain_accession",
        "full_score",
        "full_evalue",
        "domain_score",
        "independent_evalue",
        "domain_number",
        "domain_count",
        "ali_from",
        "ali_to",
        "env_from",
        "env_to",
        "accuracy",
        "source_file",
        "description",
        "3di_sequence",
    ]

    if output_path:
        output_path.parent.mkdir(parents=True, exist_ok=True)

    handle = output_path.open("w", newline="") if output_path else sys.stdout
    try:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for rank, hit in enumerate(hits, start=1):
            writer.writerow(
                [
                    rank,
                    hit.target_name,
                    hit.domain_name,
                    hit.domain_accession,
                    hit.full_score,
                    hit.full_evalue,
                    hit.domain_score,
                    hit.independent_evalue,
                    hit.domain_number,
                    hit.domain_count,
                    hit.ali_from,
                    hit.ali_to,
                    hit.env_from,
                    hit.env_to,
                    hit.accuracy,
                    hit.source_file,
                    hit.description,
                    hit.three_di_sequence,
                ]
            )
    finally:
        if output_path:
            handle.close()


def main() -> int:
    args = parse_args()
    if args.n < 1:
        raise SystemExit("n must be at least 1")

    paths = domtblout_files(args.input, args.include_plain)
    if not paths:
        raise SystemExit(f"No .with_3di.domtblout files found in {args.input}")

    matching_hits = [
        hit for hit in iter_hits(paths) if domain_matches(hit, args.domain)
    ]
    if args.unique_target:
        matching_hits = best_per_target(matching_hits, args.score)

    top_hits = select_top_hits(matching_hits, args.score, args.n)
    if not top_hits:
        searched_files = ", ".join(path.name for path in paths)
        raise SystemExit(
            f"No hits found for domain '{args.domain}' in {args.input}. "
            f"Searched {len(paths)} file(s): {searched_files}"
        )

    write_hits(top_hits, args.output)
    if args.output:
        print(f"Wrote {len(top_hits)} hit(s) to {args.output}")
    elif len(top_hits) < args.n:
        print(
            f"# Only found {len(top_hits)} hit(s) for {args.domain}; requested {args.n}.",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
