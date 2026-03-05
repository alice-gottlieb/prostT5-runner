# Installing ProstT5 and Foldseek on Hoffman2 #
Important: When installing packages, set `uv pip install --only-binary=:all: [package name]`. This is necessary for Hoffman to install and run the packages correctly

Important 2: Run all python files with `uv run python [name of file].py`. This is the only way I've found so far to get uv to correctly apply dependencies.

Run installs as follows:
```bash
uv pip install --only-binary=:all: numpy
uv pip install --only-binary=:all: transformers
uv pip install --only-binary=:all: sentencepiece
uv pip install --only-binary=:all: protobuf
uv pip install --only-binary=:all: pandas
```

If there are problems with python, set version to `uv python install 3.12`

1.  Install uv
2. Activate uv venv provided
3. Install foldseek 
`wget https://mmseqs.com/foldseek/foldseek-linux-gpu.tar.gz; tar xvfz foldseek-linux-gpu.tar.gz; export PATH=$(pwd)/foldseek/bin/:$PATH`

## batch_3di_foldseek.py ##

Uses foldseek for both 3Di prediction (via integrated ProstT5) and all-vs-all search, without independently loading ProstT5 from HuggingFace. Requires foldseek v9+ with ProstT5 support.

```bash
# Basic usage (auto-downloads ProstT5 weights, uses CPU)
uv run python batch_3di_foldseek.py assemblies.txt

# Custom output directory
uv run python batch_3di_foldseek.py assemblies.txt -o my_results

# With GPU acceleration
uv run python batch_3di_foldseek.py assemblies.txt --gpu

# Multi-threaded with GPU
uv run python batch_3di_foldseek.py assemblies.txt --gpu --threads 8

# With pre-downloaded ProstT5 weights (skips download on subsequent runs)
uv run python batch_3di_foldseek.py assemblies.txt --prostt5-weights /path/to/prostt5

# With a specific foldseek binary
uv run python batch_3di_foldseek.py assemblies.txt --foldseek-path /opt/foldseek/bin/foldseek

# With NCBI API key (faster downloads, 10 req/s vs 3 req/s)
uv run python batch_3di_foldseek.py assemblies.txt --ncbi-api-key YOUR_KEY

# Only generate 3Di codes, skip the all-vs-all search
uv run python batch_3di_foldseek.py assemblies.txt --skip-foldseek

# Pass extra arguments to foldseek search (e.g. e-value cutoff)
uv run python batch_3di_foldseek.py assemblies.txt --foldseek-args -e 0.001

# Full example with multiple options
uv run python batch_3di_foldseek.py assemblies.txt \
    -o results --gpu --threads 8 \
    --ncbi-api-key YOUR_KEY \
    --foldseek-args -e 0.001 --max-seqs 1000
```
