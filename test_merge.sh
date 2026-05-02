#!/usr/bin/env bash
#
# Smoke test for collect_3di_dbs.py / merge_two_runs.py.
#
# Usage:
#   ./test_merge.sh <run_a> <run_b> <output_dir> [foldseek_bin]
#
# Each <run_*> should be the OUTPUT_BASE of a finished submit_3di_array.sh run
# (i.e. contains task_*/foldseek_db/sequenceDB*). <output_dir> must not exist
# yet (or will be overwritten). Default foldseek binary is `foldseek` on PATH.

set -euo pipefail

if [ $# -lt 3 ]; then
    echo "Usage: $0 <run_a> <run_b> <output_dir> [foldseek_bin]" >&2
    exit 1
fi

RUN_A="$1"
RUN_B="$2"
OUT="$3"
FS="${4:-foldseek}"

REPO_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$REPO_DIR"

ORIG_HASHES="$(mktemp -t merge_test_orig.XXXXXX.sha256)"
MERGED_HEADERS="$(mktemp -t merge_test_merged.XXXXXX.txt)"
SOURCE_HEADERS="$(mktemp -t merge_test_source.XXXXXX.txt)"
trap 'rm -f "$ORIG_HASHES" "$MERGED_HEADERS" "$SOURCE_HEADERS"' EXIT

section() { printf "\n=== %s ===\n" "$1"; }

# 1) Discover all DBs
section "1) Discovering DBs"
python -c "
from pathlib import Path
from collect_3di_dbs import discover_completed_dirs, has_full_db
for run in ['$RUN_A', '$RUN_B']:
    dirs = discover_completed_dirs(Path(run))
    print(f'\n{run}: {len(dirs)} completed')
    for d in dirs:
        print(f'  {\"OK         \" if has_full_db(d) else \"FASTA-only \"}{d}')
"

# 2) Snapshot originals, then run the merge
section "2) Hashing originals + running merge"
find "$RUN_A" "$RUN_B" -path '*/foldseek_db/sequenceDB*' -type f \
    | sort | xargs shasum -a 256 > "$ORIG_HASHES"
echo "Hashed $(wc -l < "$ORIG_HASHES") source files"

python merge_two_runs.py "$RUN_A" "$RUN_B" "$OUT" -f "$FS"

# 3) Confirm originals were untouched
section "3) Verifying originals unmodified"
if shasum -a 256 -c "$ORIG_HASHES" | grep -v ': OK$'; then
    echo "FAIL: at least one original file changed" >&2
    exit 1
fi
echo "All originals intact ✓"

# 4) Verify merged DB content
section "4a) Index line counts"
SRC_COUNT=$(find "$RUN_A" "$RUN_B" -path '*/foldseek_db/sequenceDB.index' \
    -exec wc -l {} + | tail -n1 | awk '{print $1}')
MERGED_COUNT=$(wc -l < "$OUT/mergedDB.index")
echo "source: $SRC_COUNT   merged: $MERGED_COUNT"
if [ "$SRC_COUNT" != "$MERGED_COUNT" ]; then
    echo "FAIL: merged index count does not match source total" >&2
    exit 1
fi
echo "Counts match ✓"

section "4b) 3Di FASTA header diff"
"$FS" convert2fasta "$OUT/mergedDB_ss" "$OUT/merged_3di.fasta"
grep '^>' "$OUT/merged_3di.fasta" | sort -u > "$MERGED_HEADERS"
find "$RUN_A" "$RUN_B" -name all_sequences_3di.fasta \
    -exec grep -h '^>' {} + | sort -u > "$SOURCE_HEADERS"
if ! diff "$SOURCE_HEADERS" "$MERGED_HEADERS"; then
    echo "FAIL: header sets differ" >&2
    exit 1
fi
echo "Headers match ✓"

section "4c) 3Di sequence round-trip"
python - <<EOF
import glob
from Bio import SeqIO

src = {}
for f in glob.glob("$RUN_A/**/all_sequences_3di.fasta", recursive=True) + \
         glob.glob("$RUN_B/**/all_sequences_3di.fasta", recursive=True):
    for r in SeqIO.parse(f, "fasta"):
        src[r.id] = str(r.seq)

mismatches = 0
for r in SeqIO.parse("$OUT/merged_3di.fasta", "fasta"):
    if src.get(r.id) != str(r.seq):
        mismatches += 1
        if mismatches <= 5:
            print(f"MISMATCH {r.id}")

print(f"{mismatches} mismatches across {len(src)} sequences")
raise SystemExit(1 if mismatches else 0)
EOF

section "4d) Merged metadata.json"
python - <<EOF
import json
m = json.load(open("$OUT/metadata.json"))
print("source runs:", len(m["source_runs"]))
print("unique downloaded assemblies:", len(m["assemblies_downloaded"]))
print("total num_sequences:", m["num_sequences"])
EOF

section "All checks passed ✓"
