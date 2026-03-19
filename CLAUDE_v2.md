# CLAUDE.md — Workflow-Oriented Guide

## What This Project Does

prostT5-runner converts bacterial protein sequences (FASTA format) into 3Di structural alphabet codes using the ProstT5 transformer model, then runs all-vs-all structural similarity searches with Foldseek. Designed for HPC clusters with GPU acceleration.

**Conceptual pipeline:**
```
NCBI accession list → download FASTAs → ProstT5 3Di prediction → Foldseek search → results TSV
```

## Tech Stack

| Layer | Tool |
|-------|------|
| Language | Python 3.12+ |
| Package manager | `uv` |
| ML model | ProstT5 (via Hugging Face `transformers` or Foldseek built-in) |
| Structure search | Foldseek (GPU version) |
| Data | NCBI Datasets API, BioPython |
| HPC scheduler | SGE (Hoffman2) |
| GPU | CUDA 12 (H200, RTX, GTX tested) |

## Common Workflows

### 1. Full pipeline (recommended)

```bash
# assemblies.txt: one NCBI accession per line
uv run python batch_3di_foldseek.py assemblies.txt --gpu --threads 8
```

Output in `results/` (gitignored):
- `*.3di.fasta` — 3Di structural sequences
- `*.foldseek_results.tsv` — similarity hits

### 2. Large-scale parallel run (HPC)

```bash
# Step 1: fetch accessions and split into chunks
python fetch_bacterial_accessions.py --num-genomes 10000 --chunks 50
python split_accessions.py reference_genomes.txt --chunks 50

# Step 2: submit job arrays
qsub submit_download_genomes.sh
qsub submit_3di_array.sh

# Step 3: check progress
python check_progress.py

# Step 4: deduplicate and consolidate
python remove_duplicates_from_chunks.py
python collect_accessions.py
```

### 3. Analyze results

```bash
python analyze_results.py results/ overview
python analyze_results.py results/ lookup <protein_id>
python analyze_results.py results/ similar <protein_id> --evalue 1e-10 --output-fasta hits.fasta
```

### 4. Quick demo / sanity check

```bash
uv run python quick-start-Translation.py  # AA → 3Di → AA round-trip
uv run python quick-start-test.py          # embedding generation
```

## Environment Setup

```bash
# 1. Install uv
pip install uv

# 2. Install Python 3.12
uv python install 3.12

# 3. Install Python packages
#    On HPC: --only-binary avoids compilation failures
uv pip install --only-binary=:all: numpy transformers sentencepiece protobuf pandas torch

# 4. Install Foldseek (GPU build)
wget https://mmseqs.com/foldseek/foldseek-linux-gpu.tar.gz
tar xvfz foldseek-linux-gpu.tar.gz
export PATH=$(pwd)/foldseek/bin/:$PATH

# 5. Optional: NCBI API key for faster downloads
export NCBI_API_KEY=your_key_here
# or pass as --ncbi-api-key flag
```

## Codebase Map

```
batch_3di_foldseek.py          ← PRIMARY entry point (Foldseek-integrated ProstT5)
batch_3di_search.py            ← Alternative (direct Hugging Face ProstT5 loading)
download_and_process_protein.py ← Single-genome utility

analyze_results.py             ← Query/export search results
check_progress.py              ← Monitor parallel job completion
benchmark_foldseek.py          ← Performance profiling

fetch_bacterial_accessions.py  ← Pull RefSeq accession list from NCBI
split_accessions.py            ← Divide work for job arrays
collect_accessions.py          ← Aggregate completed accessions
remove_duplicates_from_chunks.py
chunk_remaining.py             ← Re-queue failed chunks
rechunk_deduped.py
quick_validate.py

submit_3di_array.sh            ← SGE job array for 3Di generation
submit_3di_parallel.sh         ← Single-node parallelism
submit_download_genomes.sh     ← SGE job array for downloads
flatten_genomes.sh             ← Post-download file organization

quick-start-Translation.py     ← Demo: bidirectional AA ↔ 3Di
quick-start-test.py            ← Demo: embeddings
```

## Gotchas

- `batch_3di_foldseek.py` is faster than `batch_3di_search.py` because Foldseek has ProstT5 compiled in — prefer it unless you need direct model access.
- The `foldseek/` binary directory and all `results*/` directories are gitignored.
- Store your NCBI API key in `ncbi_key.txt` (also gitignored) or pass via flag.
- Benchmark data for H200 (32GB), RTX (16GB), GTX (16GB) is in `benchmarks/`.
