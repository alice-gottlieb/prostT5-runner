# CLAUDE.md

## Project Overview

**prostT5-runner** is a bioinformatics pipeline for large-scale bacterial protein structure analysis. It downloads bacterial reference genome protein sequences from NCBI, converts amino acid sequences to 3Di structural codes using ProstT5, and performs all-vs-all structural similarity searches using Foldseek.

## Commands

### Run / Entry Points

```bash
# Primary pipeline: download + 3Di prediction + similarity search
uv run python batch_3di_foldseek.py assemblies.txt

# With GPU and threads
uv run python batch_3di_foldseek.py assemblies.txt --gpu --threads 8

# Skip the Foldseek search step (only generate 3Di codes)
uv run python batch_3di_foldseek.py assemblies.txt --skip-foldseek

# Alternative: use ProstT5 directly via Hugging Face (slower)
uv run python batch_3di_search.py assemblies.txt --gpu

# Quick demos
uv run python quick-start-Translation.py
uv run python quick-start-test.py
```

### Analysis

```bash
uv run python analyze_results.py results_dir overview
uv run python analyze_results.py results_dir lookup WP_004995605.1
uv run python analyze_results.py results_dir similar WP_004995605.1 --evalue 1e-10
uv run python analyze_results.py results_dir similar WP_004995605.1 --output-fasta hits.fasta
```

### HPC (Hoffman2 / SGE)

```bash
qsub submit_3di_array.sh         # Distributed 3Di generation
qsub submit_download_genomes.sh  # Distributed genome downloading
```

### Benchmarking

```bash
uv run python benchmark_foldseek.py assemblies.txt \
    --foldseek-path foldseek/bin/foldseek \
    --prostt5-weights /path/to/prostt5 \
    --gpu --thread-counts 1 2 4 8
```

### Dependencies

```bash
# Install Python deps (use --only-binary on HPC to avoid compilation)
uv pip install --only-binary=:all: numpy transformers sentencepiece protobuf pandas

# Install Foldseek (GPU version)
wget https://mmseqs.com/foldseek/foldseek-linux-gpu.tar.gz
tar xvfz foldseek-linux-gpu.tar.gz
export PATH=$(pwd)/foldseek/bin/:$PATH
```

## Architecture

- **`batch_3di_foldseek.py`** — Main pipeline using Foldseek's built-in ProstT5 integration (preferred)
- **`batch_3di_search.py`** — Alternative pipeline loading ProstT5 directly via Hugging Face
- **`download_and_process_protein.py`** — Single-genome processing
- **`analyze_results.py`** — Query and export similarity search results
- **`check_progress.py`** — Monitor progress across parallel job chunks
- **`fetch_bacterial_accessions.py`** / **`split_accessions.py`** — Data management for parallel runs

## Key Notes

- Python 3.12+ required; use `uv` as the package manager
- GPU support via CUDA 12 (tested on H200, RTX, GTX)
- On Hoffman2, always use `--only-binary=:all:` to avoid compilation failures
- NCBI API key (`--ncbi-api-key`) raises download rate limit from 3 to 10 req/s
- Results and logs are gitignored; `foldseek/` binary is also excluded
