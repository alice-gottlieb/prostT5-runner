#!/usr/bin/env python3
"""
Download Pfam release files from EBI FTP for downstream HMMER scanning of
bacterial reference genomes.

Pfam moved from pfam.xfam.org to EBI; current files live at:
    https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/

Usage:
    python download_pfam.py --output-dir pfam_data
    python download_pfam.py --output-dir pfam_data --files hmm,dat --verify
"""

import argparse
import hashlib
import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path

import requests


PFAM_BASE = "https://ftp.ebi.ac.uk/pub/databases/Pfam"

# Logical name -> remote filename. The "always" group is small and always fetched.
FILES = {
    "hmm":   "Pfam-A.hmm.gz",
    "dat":   "Pfam-A.hmm.dat.gz",
    "full":  "Pfam-A.full.gz",
    "fasta": "Pfam-A.fasta.gz",
}
ALWAYS_FILES = ["relnotes.txt", "Pfam.version.gz", "md5_checksums"]

CHUNK_SIZE = 1 << 20  # 1 MiB
PROGRESS_EVERY_MB = 100


def release_url(release: str) -> str:
    return f"{PFAM_BASE}/{release}"


def download_file(url: str, dest: Path, skip_existing: bool) -> Path:
    """Stream a single file to disk. Skips if already present and skip_existing."""
    if skip_existing and dest.exists() and dest.stat().st_size > 0:
        print(f"  skip (exists): {dest.name}")
        return dest

    print(f"  GET {url}")
    start = time.time()
    with requests.get(url, stream=True, timeout=60) as resp:
        resp.raise_for_status()
        total = int(resp.headers.get("Content-Length", 0))
        bytes_done = 0
        next_report = PROGRESS_EVERY_MB * (1 << 20)

        tmp = dest.with_suffix(dest.suffix + ".part")
        with open(tmp, "wb") as f:
            for chunk in resp.iter_content(chunk_size=CHUNK_SIZE):
                if not chunk:
                    continue
                f.write(chunk)
                bytes_done += len(chunk)
                if bytes_done >= next_report:
                    mb = bytes_done / (1 << 20)
                    pct = f" ({100 * bytes_done / total:.1f}%)" if total else ""
                    print(f"    {mb:.0f} MB{pct}")
                    next_report += PROGRESS_EVERY_MB * (1 << 20)

        tmp.rename(dest)

    elapsed = time.time() - start
    mb = dest.stat().st_size / (1 << 20)
    rate = mb / elapsed if elapsed > 0 else 0
    print(f"  done: {dest.name} ({mb:.1f} MB in {elapsed:.0f}s, {rate:.1f} MB/s)")
    return dest


def md5sum(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(CHUNK_SIZE), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_md5_checksums(path: Path) -> dict:
    """Pfam's md5_checksums file: one '<hash>  <filename>' per line."""
    sums = {}
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2:
                sums[parts[1]] = parts[0]
    return sums


def verify(dest: Path, expected_md5: str) -> bool:
    actual = md5sum(dest)
    ok = actual == expected_md5
    status = "OK" if ok else "MISMATCH"
    print(f"  md5 {status}: {dest.name} (expected {expected_md5}, got {actual})")
    return ok


def read_release_version(notes_path: Path) -> str:
    """Pull the first 'Release X.Y' line out of relnotes.txt."""
    if not notes_path.exists():
        return "unknown"
    with open(notes_path) as f:
        for line in f:
            line = line.strip()
            if line.lower().startswith("release"):
                return line
    return "unknown"


def parse_args():
    p = argparse.ArgumentParser(description="Download Pfam release files from EBI.")
    p.add_argument("--output-dir", default="pfam_data",
                   help="Directory to write Pfam files into (default: pfam_data)")
    p.add_argument("--release", default="current_release",
                   help="Release subdirectory on the FTP (default: current_release)")
    p.add_argument("--files", default="hmm,dat,full,fasta",
                   help=f"Comma-separated subset of {sorted(FILES)} (default: all)")
    p.add_argument("--skip-existing", action="store_true",
                   help="Skip files that already exist in --output-dir")
    p.add_argument("--verify", action="store_true",
                   help="Verify md5 of each downloaded file against Pfam's md5_checksums")
    return p.parse_args()


def main():
    args = parse_args()

    requested = [s.strip() for s in args.files.split(",") if s.strip()]
    unknown = [s for s in requested if s not in FILES]
    if unknown:
        raise SystemExit(f"Unknown --files entries: {unknown}. Valid: {sorted(FILES)}")

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    base = release_url(args.release)

    print(f"Pfam release URL: {base}")
    print(f"Output directory: {out_dir.resolve()}")
    print(f"Files: {requested + ALWAYS_FILES}")

    # Fetch metadata + checksums first.
    for name in ALWAYS_FILES:
        download_file(f"{base}/{name}", out_dir / name, args.skip_existing)

    # Fetch the bulk files.
    fetched = []
    for key in requested:
        fname = FILES[key]
        dest = download_file(f"{base}/{fname}", out_dir / fname, args.skip_existing)
        fetched.append(fname)

    # Optional checksum verification.
    if args.verify:
        sums = parse_md5_checksums(out_dir / "md5_checksums")
        all_ok = True
        for fname in fetched:
            expected = sums.get(fname)
            if expected is None:
                print(f"  md5 MISSING: no entry for {fname} in md5_checksums")
                all_ok = False
                continue
            if not verify(out_dir / fname, expected):
                all_ok = False
        if not all_ok:
            raise SystemExit("md5 verification failed for one or more files")

    # Record what we did.
    metadata = {
        "release_url": base,
        "release_version": read_release_version(out_dir / "relnotes.txt"),
        "downloaded_at_utc": datetime.now(timezone.utc).isoformat(),
        "files": [
            {"name": f, "bytes": (out_dir / f).stat().st_size}
            for f in fetched + ALWAYS_FILES
            if (out_dir / f).exists()
        ],
        "verified_md5": bool(args.verify),
    }
    with open(out_dir / "download_metadata.json", "w") as f:
        json.dump(metadata, f, indent=2)

    print(f"\nDone. Wrote metadata to {out_dir / 'download_metadata.json'}")


if __name__ == "__main__":
    main()
