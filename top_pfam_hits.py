#!/usr/bin/env python3
"""Find the top HMMER hits for a Pfam domain in pfam_hits2-1by1.

Examples:
    uv run python top_pfam_hits.py AAA 10
    uv run python top_pfam_hits.py PF00004 10 --unique-target
    uv run python top_pfam_hits.py PF00004.36 10 --score full
    uv run python top_pfam_hits.py PF00198 5 -i pfam_hits2-1by1/GCF_000005845.2.with_3di.domtblout
    uv run python top_pfam_hits.py --all-domains 10 -o all_top_hits.tsv
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


@dataclass(frozen=True)
class PfamDomain:
    accession: str
    pfam_id: str
    description: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Find the top N rows for a Pfam domain across HMMER domtblout files "
            "in pfam_hits2-1by1."
        )
    )
    parser.add_argument(
        "domain",
        nargs="?",
        help=(
            "Pfam domain query name or accession. Examples: AAA, PF00004, "
            "PF00004.36."
        ),
    )
    parser.add_argument("n", nargs="?", type=int, help="Number of top hits to print.")
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
    parser.add_argument(
        "--all-domains",
        action="store_true",
        help=(
            "Find the top N hits for every Pfam domain listed in "
            "--pfam-domains."
        ),
    )
    parser.add_argument(
        "--pfam-domains",
        type=Path,
        default=Path("pfam_domains.tsv"),
        help=(
            "TSV containing pfam_accession, pfam_id, and description columns. "
            "Used with --all-domains."
        ),
    )
    parser.add_argument(
        "--progress-interval",
        type=int,
        default=1000,
        help=(
            "Print and flush progress every N parsed hit rows and every N Pfam "
            "domains in --all-domains mode. Use 0 to disable periodic updates."
        ),
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


def log_progress(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


def iter_hits(paths: Iterable[Path], progress_interval: int = 1000) -> Iterable[PfamHit]:
    for path in paths:
        log_progress(f"Reading hits from {path}")
        parsed_rows = 0
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
                        flush=True,
                    )
                    continue
                if hit is not None:
                    parsed_rows += 1
                    if progress_interval > 0 and parsed_rows % progress_interval == 0:
                        log_progress(f"  parsed {parsed_rows} hit row(s) from {path}")
                    yield hit
        log_progress(f"Finished {path}: parsed {parsed_rows} hit row(s)")


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


def output_columns(prefix_domain: bool = False) -> list[str]:
    columns = []
    if prefix_domain:
        columns.extend(
            [
                "query_pfam_accession",
                "query_pfam_id",
                "query_pfam_description",
            ]
        )
    columns.extend(
        [
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
    )
    return columns


def hit_values(hit: PfamHit, rank: int) -> list[object]:
    return [
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


def open_output(output_path: Path | None):
    if output_path:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        return output_path.open("w", newline="")
    return sys.stdout


def write_hits(hits: list[PfamHit], output_path: Path | None) -> None:
    handle = open_output(output_path)
    try:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(output_columns())
        for rank, hit in enumerate(hits, start=1):
            writer.writerow(hit_values(hit, rank))
    finally:
        if output_path:
            handle.close()


def read_pfam_domains(path: Path) -> list[PfamDomain]:
    if not path.is_file():
        raise SystemExit(f"Pfam domains file not found: {path}")

    domains: list[PfamDomain] = []
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row_number, row in enumerate(reader, start=1):
            if not row:
                continue
            if row_number == 1 and row[0] == "pfam_accession":
                continue
            if len(row) < 3:
                print(
                    f"WARN: skipped malformed Pfam domain row {path}:{row_number}",
                    file=sys.stderr,
                    flush=True,
                )
                continue
            domains.append(PfamDomain(row[0], row[1], row[2]))

    return domains


def domain_key_lookup(domains: Iterable[PfamDomain]) -> dict[str, str]:
    lookup: dict[str, str] = {}
    for domain in domains:
        lookup[normalize_domain(domain.accession)] = domain.accession
        lookup[normalize_domain(domain.pfam_id)] = domain.accession
    return lookup


def hit_domain_accession(
    hit: PfamHit,
    lookup: dict[str, str],
) -> str | None:
    return lookup.get(normalize_domain(hit.domain_accession)) or lookup.get(
        normalize_domain(hit.domain_name)
    )


def trim_candidates(
    candidates: list[PfamHit],
    score: str,
    n: int,
) -> list[PfamHit]:
    return select_top_hits(candidates, score, n)


def collect_all_domain_candidates(
    hits: Iterable[PfamHit],
    domains: list[PfamDomain],
    score: str,
    n: int,
    unique_target: bool,
) -> dict[str, list[PfamHit]] | dict[str, dict[str, PfamHit]]:
    lookup = domain_key_lookup(domains)

    if unique_target:
        unique_candidates: dict[str, dict[str, PfamHit]] = {}
        for hit in hits:
            accession = hit_domain_accession(hit, lookup)
            if accession is None:
                continue
            target_hits = unique_candidates.setdefault(accession, {})
            current = target_hits.get(hit.target_name)
            if current is None or hit_sort_key(hit, score) < hit_sort_key(current, score):
                target_hits[hit.target_name] = hit
        return unique_candidates

    candidates: dict[str, list[PfamHit]] = {}
    trim_at = max(n * 4, 100)
    keep_at = max(n * 2, n)
    for hit in hits:
        accession = hit_domain_accession(hit, lookup)
        if accession is None:
            continue
        domain_candidates = candidates.setdefault(accession, [])
        domain_candidates.append(hit)
        if len(domain_candidates) > trim_at:
            candidates[accession] = trim_candidates(domain_candidates, score, keep_at)
    return candidates


def write_all_domain_hits(
    domains: list[PfamDomain],
    candidates: dict[str, list[PfamHit]] | dict[str, dict[str, PfamHit]],
    score: str,
    n: int,
    unique_target: bool,
    output_path: Path | None,
    progress_interval: int,
) -> tuple[int, int]:
    domains_with_hits = 0
    rows_written = 0
    handle = open_output(output_path)
    try:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(output_columns(prefix_domain=True))
        handle.flush()
        for domain_number, domain in enumerate(domains, start=1):
            if progress_interval > 0 and domain_number % progress_interval == 0:
                log_progress(
                    f"Processed {domain_number}/{len(domains)} Pfam domain(s); "
                    f"wrote {rows_written} hit row(s)"
                )
            if unique_target:
                domain_candidates = candidates.get(domain.accession, {})
                hits = list(domain_candidates.values())
            else:
                domain_candidates = candidates.get(domain.accession, [])
                hits = domain_candidates
            top_hits = select_top_hits(hits, score, n)
            if not top_hits:
                handle.flush()
                continue

            domains_with_hits += 1
            rows_written += len(top_hits)
            for rank, hit in enumerate(top_hits, start=1):
                writer.writerow(
                    [
                        domain.accession,
                        domain.pfam_id,
                        domain.description,
                        *hit_values(hit, rank),
                    ]
                )
            handle.flush()
    finally:
        if output_path:
            handle.close()

    return domains_with_hits, rows_written


def parse_n(args: argparse.Namespace) -> int:
    if args.all_domains and args.n is None and args.domain is not None:
        try:
            args.n = int(args.domain)
        except ValueError:
            pass
        else:
            args.domain = None

    if args.n is None:
        raise SystemExit("n is required")
    if args.n < 1:
        raise SystemExit("n must be at least 1")
    if args.progress_interval < 0:
        raise SystemExit("--progress-interval must be >= 0")
    return args.n


def main() -> int:
    args = parse_args()
    n = parse_n(args)
    if not args.all_domains and args.domain is None:
        raise SystemExit("domain is required unless --all-domains is used")

    paths = domtblout_files(args.input, args.include_plain)
    if not paths:
        raise SystemExit(f"No .with_3di.domtblout files found in {args.input}")

    if args.all_domains:
        log_progress(f"Reading Pfam domain list from {args.pfam_domains}")
        domains = read_pfam_domains(args.pfam_domains)
        log_progress(f"Loaded {len(domains)} Pfam domain(s)")
        log_progress(
            f"Collecting top hit candidates from {len(paths)} domtblout file(s)"
        )
        candidates = collect_all_domain_candidates(
            hits=iter_hits(paths, progress_interval=args.progress_interval),
            domains=domains,
            score=args.score,
            n=n,
            unique_target=args.unique_target,
        )
        log_progress(f"Collected candidates for {len(candidates)} Pfam domain(s)")
        log_progress("Writing top hits for each Pfam domain")
        domains_with_hits, rows_written = write_all_domain_hits(
            domains=domains,
            candidates=candidates,
            score=args.score,
            n=n,
            unique_target=args.unique_target,
            output_path=args.output,
            progress_interval=args.progress_interval,
        )
        log_progress(
            f"Wrote {rows_written} hit row(s) for {domains_with_hits}/"
            f"{len(domains)} Pfam domain(s)."
        )
        return 0

    log_progress(f"Searching for top {n} hit(s) for {args.domain}")
    matching_hits = [
        hit
        for hit in iter_hits(paths, progress_interval=args.progress_interval)
        if domain_matches(hit, args.domain)
    ]
    log_progress(f"Found {len(matching_hits)} matching hit row(s) for {args.domain}")
    if args.unique_target:
        matching_hits = best_per_target(matching_hits, args.score)
        log_progress(f"Kept {len(matching_hits)} best unique target hit row(s)")

    top_hits = select_top_hits(matching_hits, args.score, n)
    if not top_hits:
        searched_files = ", ".join(path.name for path in paths)
        raise SystemExit(
            f"No hits found for domain '{args.domain}' in {args.input}. "
            f"Searched {len(paths)} file(s): {searched_files}"
        )

    write_hits(top_hits, args.output)
    if args.output:
        print(f"Wrote {len(top_hits)} hit(s) to {args.output}", flush=True)
    elif len(top_hits) < n:
        print(
            f"# Only found {len(top_hits)} hit(s) for {args.domain}; requested {n}.",
            file=sys.stderr,
            flush=True,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
