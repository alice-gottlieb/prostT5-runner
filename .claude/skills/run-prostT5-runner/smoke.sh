#!/usr/bin/env bash
# Smoke driver for the foldseek_topn_pfam pipeline (the actively-developed
# core of prostT5-runner). Two phases:
#   1. Run the repo's output-integrity test suite (logic groups A-D run
#      anywhere; group E shells out to the real foldseek binary).
#   2. Run the real CLI end-to-end on the checked-in sample input + target DB
#      and assert the two output files exist and are non-empty.
#
# Run from the repo root:
#   bash .claude/skills/run-prostT5-runner/smoke.sh
#
# Exits 0 only if the test suite passes AND the CLI produces hits.tsv +
# genome_counts.tsv. Honors $FOLDSEEK (defaults to ~/foldseek/bin/foldseek).
set -euo pipefail

REPO="$(git rev-parse --show-toplevel)"
cd "$REPO"

FOLDSEEK="${FOLDSEEK:-$HOME/foldseek/bin/foldseek}"
TARGET_DB="results_fs_only_full_fs_compare/foldseek_db/sequenceDB"
OUT="tmp/fstopn_smoke"   # repo-local scratch dir

echo "=== prostT5-runner smoke driver ==="
echo "repo:     $REPO"
echo "foldseek: $FOLDSEEK"
[ -x "$FOLDSEEK" ] || { echo "FATAL: foldseek not found/executable at $FOLDSEEK"; exit 1; }
[ -f "${TARGET_DB}.dbtype" ] || { echo "FATAL: target DB missing at $TARGET_DB"; exit 1; }

echo
echo "--- Phase 1: output-integrity test suite ---"
FOLDSEEK="$FOLDSEEK" uv run python test_foldseek_topn_pfam.py | tail -1

echo
echo "--- Phase 2: real CLI end-to-end ---"
mkdir -p tmp
rm -rf "$OUT"
uv run python foldseek_topn_pfam.py \
    pf00198_top5.tsv "$TARGET_DB" "$OUT" \
    --foldseek "$FOLDSEEK" --gpu 0 --threads 2 --max-seqs 5 \
    > "$OUT.log" 2>&1 || { echo "FATAL: CLI run failed"; tail -20 "$OUT.log"; exit 1; }

for f in hits.tsv genome_counts.tsv; do
    [ -s "$OUT/$f" ] || { echo "FATAL: $OUT/$f missing or empty"; exit 1; }
    echo "  ok: $OUT/$f ($(wc -l < "$OUT/$f" | tr -d ' ') lines)"
done

echo
echo "=== SMOKE OK ==="
