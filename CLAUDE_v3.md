# CLAUDE.md — Minimal Reference Card

## Purpose

Bacterial protein structure pipeline: NCBI accessions → ProstT5 3Di codes → Foldseek similarity search.

## Quick Commands

```bash
# Full pipeline
uv run python batch_3di_foldseek.py assemblies.txt --gpu --threads 8

# 3Di only (no search)
uv run python batch_3di_foldseek.py assemblies.txt --skip-foldseek

# Analyze results
uv run python analyze_results.py results/ overview
uv run python analyze_results.py results/ similar <protein_id> --evalue 1e-10

# HPC submission (SGE)
qsub submit_3di_array.sh

# Demos
uv run python quick-start-Translation.py
```

## Setup

```bash
uv python install 3.12
uv pip install --only-binary=:all: numpy transformers sentencepiece protobuf pandas torch
# Install Foldseek GPU binary and add to PATH
```

## Key Files

| File | Role |
|------|------|
| `batch_3di_foldseek.py` | Main entry point |
| `batch_3di_search.py` | Alt: direct HuggingFace ProstT5 |
| `analyze_results.py` | Query results |
| `check_progress.py` | Monitor HPC jobs |
| `fetch_bacterial_accessions.py` | Get NCBI accessions |
| `split_accessions.py` | Chunk for parallelism |
| `submit_3di_array.sh` | SGE job array |
| `benchmark_foldseek.py` | Performance testing |

## Notes

- Python 3.12+, `uv` package manager
- CUDA 12 GPU required for `--gpu` flag
- Use `--only-binary=:all:` on HPC (Hoffman2) to avoid build failures
- `results*/`, `logs/`, `foldseek/`, `ncbi_key.txt` are all gitignored
