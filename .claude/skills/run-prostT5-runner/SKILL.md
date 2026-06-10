---
name: run-prostT5-runner
description: Build, run, test, and smoke-check the prostT5-runner foldseek/3Di pipelines. Use when asked to run, drive, smoke-test, or verify the foldseek_topn_pfam pipeline or batch_3di_foldseek 3Di-prediction tool, or to confirm a change to these scripts works end-to-end.
---

# Run prostT5-runner

prostT5-runner is a collection of bioinformatics CLI scripts that predict 3Di
structural codes for protein sequences (via Foldseek's integrated ProstT5) and
run Foldseek searches. The actively-developed core is **`foldseek_topn_pfam.py`**:
it takes a top-N Pfam-hits TSV of 3Di codes, builds a Foldseek query DB,
searches it against a target 3Di DB, and aggregates hits by genome.

The agent path is a smoke driver that runs the repo's output-integrity test
suite **and** the real CLI end-to-end, asserting the outputs exist:

```bash
bash .claude/skills/run-prostT5-runner/smoke.sh
```

All paths below are relative to the repo root.

## Prerequisites

- **uv** (package manager). Run every Python file with `uv run python …` — uv
  applies the locked deps; bare `python` will not have torch/transformers.
- **foldseek v9+** with ProstT5 support, at `~/foldseek/bin/foldseek` (a GPU
  build is fine; the smoke run uses `--gpu 0`). Override with `$FOLDSEEK`.
  Install (Linux): `wget https://mmseqs.com/foldseek/foldseek-linux-gpu.tar.gz && tar xvfz foldseek-linux-gpu.tar.gz`
- No build step — these are scripts. `uv run` resolves deps on first invocation.

## Run (agent path)

```bash
bash .claude/skills/run-prostT5-runner/smoke.sh
```

Expected tail:

```
--- Phase 1: output-integrity test suite ---
ALL TESTS PASSED

--- Phase 2: real CLI end-to-end ---
  ok: tmp/fstopn_smoke/hits.tsv (11 lines)
  ok: tmp/fstopn_smoke/genome_counts.tsv (2 lines)

=== SMOKE OK ===
```

Scratch output lands in the repo-local **`tmp/fstopn_smoke/`** (gitignored).
`hits.tsv` is one Foldseek hit per row plus query metadata + target genome;
`genome_counts.tsv` is a genome × query_id count matrix.

## Direct invocation

Run the CLI yourself with the checked-in sample input + target DB:

```bash
uv run python foldseek_topn_pfam.py \
    pf00198_top5.tsv \
    results_fs_only_full_fs_compare/foldseek_db/sequenceDB \
    tmp/fstopn_smoke \
    --foldseek ~/foldseek/bin/foldseek --gpu 0 --threads 2 --max-seqs 5
```

Larger checked-in inputs: `all_pfam_top10_hits.tsv` (add `--test` to keep only
the top 5 rows of the first 3 Pfam accessions). `uv run python foldseek_topn_pfam.py --help`
for all flags.

## Test

The output-integrity suite (proves every output cell traces back to the right
pfam/genome/column — catches mis-joins and byte-offset bugs in the hand-built
query DB):

```bash
uv run python test_foldseek_topn_pfam.py
```

Groups A–D are pure data-logic and run anywhere; group E is end-to-end and
auto-skips if foldseek or the target DB is absent. With foldseek at
`~/foldseek/bin/foldseek` and the checked-in target DB present, group E runs
(0 skips). The merge pipeline has a separate test:
`uv run python test_merge.py <run_a> <run_b> <out>` (needs two finished array runs).

## The other entry point: batch_3di_foldseek.py

`batch_3di_foldseek.py` (the README's headline tool) downloads whole-genome
`.faa` files from NCBI Datasets, converts them to 3Di with foldseek's
integrated ProstT5, then runs an all-vs-all search:

```bash
uv run python batch_3di_foldseek.py assemblies.txt -o tmp/batch_out --foldseek-path ~/foldseek/bin/foldseek
```

Not exercised by the smoke driver: it requires live NCBI network access and a
~1GB ProstT5 weights download on first run (`foldseek databases ProstT5`). Use
it when you specifically need the download→predict→search flow; otherwise the
`foldseek_topn_pfam.py` path above is the fast, fully-local way to drive the
search core.

## Gotchas

- **`pytest` is not installed** — there is no `pytest`. The tests are
  self-contained runners invoked with `uv run python test_*.py`; they print
  `ALL TESTS PASSED` / `ALL TESTS PASSED` and exit non-zero on first failure.
- **`foldseek version` prints only a git hash** (e.g. `718d4217…`), not a
  semver. That's the real binary; don't mistake it for an error.
- **`--gpu 1` is the CLI default** for `foldseek_topn_pfam.py`. On a machine
  without a usable GPU, pass `--gpu 0` (the smoke driver does).
- **Write scratch under the repo's `tmp/`**, not system `/tmp` — `tmp/**` is
  gitignored, keeping `git status` clean.
- **target_genome shows `NA`** in the sample run: the checked-in target DB IDs
  (`NP_…`, `WP_…`) don't carry a parseable genome accession, so genomes
  aggregate under `NA`. Expected for the sample DB, not a bug.

## Troubleshooting

- `error: Failed to spawn: pytest` → you tried `uv run pytest`. There is no
  pytest; run `uv run python test_foldseek_topn_pfam.py` instead.
- `FATAL: foldseek not found` from the driver → set `FOLDSEEK=/path/to/foldseek`
  before running, or install to `~/foldseek/bin/foldseek`.
- `FATAL: target DB missing` → the checked-in DB at
  `results_fs_only_full_fs_compare/foldseek_db/sequenceDB` is gitignored
  (`results**`); it must be present locally for phase 2 / group E to run.
