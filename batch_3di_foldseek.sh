#!/usr/bin/env bash
# Batch download whole-species protein FASTA (.faa) files from NCBI Datasets,
# convert all sequences to 3Di codes using foldseek's integrated ProstT5, then
# run an all-vs-all foldseek search across all sequences.
#
# This is an alternative to batch_3di_search.py that uses foldseek for both
# 3Di prediction and search, rather than independently loading ProstT5 from
# HuggingFace via Python.
#
# Input:  a text file with one NCBI genome assembly accession (GCF_* or GCA_*)
#         per line
# Output: per-assembly .faa files, a combined 3Di FASTA, and foldseek results
#
# Requires: foldseek (v9+ with ProstT5 support), curl, unzip, jq (optional, for metadata)
#
# Example usages:
#
#     # Basic usage (auto-downloads ProstT5 weights, uses CPU)
#     bash batch_3di_foldseek.sh assemblies.txt
#
#     # Custom output directory
#     bash batch_3di_foldseek.sh assemblies.txt -o my_results
#
#     # With GPU acceleration
#     bash batch_3di_foldseek.sh assemblies.txt --gpu
#
#     # Multi-threaded with GPU
#     bash batch_3di_foldseek.sh assemblies.txt --gpu --threads 8
#
#     # With pre-downloaded ProstT5 weights
#     bash batch_3di_foldseek.sh assemblies.txt --prostt5-weights /path/to/prostt5
#
#     # With a specific foldseek binary
#     bash batch_3di_foldseek.sh assemblies.txt --foldseek-path /opt/foldseek/bin/foldseek
#
#     # With NCBI API key (faster downloads, 10 req/s vs 3 req/s)
#     bash batch_3di_foldseek.sh assemblies.txt --ncbi-api-key YOUR_KEY
#
#     # Only generate 3Di codes, skip the all-vs-all search
#     bash batch_3di_foldseek.sh assemblies.txt --skip-foldseek
#
#     # Pass extra arguments to foldseek search (e.g. e-value cutoff)
#     bash batch_3di_foldseek.sh assemblies.txt --foldseek-args -e 0.001
#
#     # Full example with multiple options
#     bash batch_3di_foldseek.sh assemblies.txt \
#         -o results --gpu --threads 8 \
#         --ncbi-api-key YOUR_KEY \
#         --foldseek-args -e 0.001 --max-seqs 1000

set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
OUTPUT_DIR="batch_output"
FASTA_DIR=""
NCBI_API_KEY=""
FOLDSEEK_PATH=""
PROSTT5_WEIGHTS=""
USE_GPU=false
PROSTT5_SPLIT_LENGTH=""
THREADS=1
SKIP_FOLDSEEK=false
FOLDSEEK_ARGS=()
ACCESSION_FILE=""

# ---------------------------------------------------------------------------
# Usage
# ---------------------------------------------------------------------------
usage() {
    cat <<'EOF'
Usage: batch_3di_foldseek.sh ACCESSION_FILE [OPTIONS]

Batch download protein FASTAs from NCBI, convert to 3Di with foldseek, and
run all-vs-all search.

Positional arguments:
  ACCESSION_FILE        Text file with one NCBI assembly accession per line

Options:
  -o, --output-dir DIR          Output directory (default: batch_output)
  --fasta-dir DIR               Directory with pre-downloaded .faa files (skips download)
  --ncbi-api-key KEY            NCBI API key for higher rate limits
  --foldseek-path PATH          Path to foldseek binary or installation directory
  --prostt5-weights PATH        Path prefix for pre-downloaded ProstT5 weights
  --gpu                         Use GPU for ProstT5 inference
  --prostt5-split-length N      Max sequence length before splitting runs
  -t, --threads N               Number of threads for foldseek (default: 1)
  --skip-foldseek               Stop after generating 3Di codes
  --foldseek-args ARGS...       Extra arguments forwarded to foldseek search
  -h, --help                    Show this help message
EOF
    exit 0
}

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            usage
            ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"; shift 2
            ;;
        --fasta-dir)
            FASTA_DIR="$2"; shift 2
            ;;
        --ncbi-api-key)
            NCBI_API_KEY="$2"; shift 2
            ;;
        --foldseek-path)
            FOLDSEEK_PATH="$2"; shift 2
            ;;
        --prostt5-weights)
            PROSTT5_WEIGHTS="$2"; shift 2
            ;;
        --gpu)
            USE_GPU=true; shift
            ;;
        --prostt5-split-length)
            PROSTT5_SPLIT_LENGTH="$2"; shift 2
            ;;
        -t|--threads)
            THREADS="$2"; shift 2
            ;;
        --skip-foldseek)
            SKIP_FOLDSEEK=true; shift
            ;;
        --foldseek-args)
            shift
            FOLDSEEK_ARGS=("$@")
            break
            ;;
        -*)
            echo "Unknown option: $1" >&2; exit 1
            ;;
        *)
            if [[ -z "$ACCESSION_FILE" ]]; then
                ACCESSION_FILE="$1"; shift
            else
                echo "Unexpected argument: $1" >&2; exit 1
            fi
            ;;
    esac
done

if [[ -z "$ACCESSION_FILE" ]]; then
    echo "Error: ACCESSION_FILE is required" >&2
    usage
fi

if [[ ! -f "$ACCESSION_FILE" ]]; then
    echo "Error: accession file not found: $ACCESSION_FILE" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Resolve foldseek binary
# ---------------------------------------------------------------------------
resolve_foldseek_bin() {
    if [[ -n "$FOLDSEEK_PATH" ]]; then
        if [[ -d "$FOLDSEEK_PATH" ]]; then
            if [[ -x "$FOLDSEEK_PATH/bin/foldseek" ]]; then
                FOLDSEEK_BIN="$FOLDSEEK_PATH/bin/foldseek"
            elif [[ -x "$FOLDSEEK_PATH/foldseek" ]]; then
                FOLDSEEK_BIN="$FOLDSEEK_PATH/foldseek"
            else
                echo "Error: foldseek binary not found in $FOLDSEEK_PATH" >&2
                exit 1
            fi
        else
            if [[ ! -x "$FOLDSEEK_PATH" ]]; then
                echo "Error: foldseek binary not found at $FOLDSEEK_PATH" >&2
                exit 1
            fi
            FOLDSEEK_BIN="$FOLDSEEK_PATH"
        fi
    else
        FOLDSEEK_BIN="foldseek"
    fi
}

# ---------------------------------------------------------------------------
# Helper: run a command with logging
# ---------------------------------------------------------------------------
run_cmd() {
    local description="$1"
    shift
    echo ""
    echo "$description: $*"
    "$@" 2>&1 | head -c 4096 || {
        echo "Command failed: $*" >&2
        return 1
    }
}

# ---------------------------------------------------------------------------
# Read accessions
# ---------------------------------------------------------------------------
read_accessions() {
    local file="$1"
    ACCESSIONS=()
    while IFS= read -r line; do
        # Trim whitespace
        line="${line#"${line%%[![:space:]]*}"}"
        line="${line%"${line##*[![:space:]]}"}"
        # Skip blank lines and comments
        [[ -z "$line" || "$line" == \#* ]] && continue
        ACCESSIONS+=("$line")
    done < "$file"
    echo "Loaded ${#ACCESSIONS[@]} assembly accessions from $file"
}

# ---------------------------------------------------------------------------
# Download a single species .faa from NCBI Datasets API
# ---------------------------------------------------------------------------
download_species_faa() {
    local acc="$1"
    local output_path="$2"
    local tmp_zip
    tmp_zip=$(mktemp "${TMPDIR:-/tmp}/ncbi_download_XXXXXX.zip")

    local url="https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/${acc}/download?include_annotation_type=PROT_FASTA"
    local curl_args=(-sS -L --fail -o "$tmp_zip" --max-time 120 -H "Accept: application/zip")
    if [[ -n "$NCBI_API_KEY" ]]; then
        curl_args+=(-H "api-key: $NCBI_API_KEY")
    fi

    if ! curl "${curl_args[@]}" "$url"; then
        echo "  [WARN] $acc: download failed" >&2
        rm -f "$tmp_zip"
        return 1
    fi

    # Extract protein.faa from the ZIP
    local faa_entry
    faa_entry=$(unzip -l "$tmp_zip" 2>/dev/null | grep 'protein\.faa' | awk '{print $NF}' | head -1)
    if [[ -z "$faa_entry" ]]; then
        echo "  [WARN] $acc: no protein.faa found in downloaded zip" >&2
        rm -f "$tmp_zip"
        return 1
    fi

    unzip -p "$tmp_zip" "$faa_entry" > "$output_path"
    rm -f "$tmp_zip"
    return 0
}

# ---------------------------------------------------------------------------
# Batch download
# ---------------------------------------------------------------------------
batch_download() {
    local fasta_dir="$1"
    mkdir -p "$fasta_dir"

    local delay=0.34
    if [[ -n "$NCBI_API_KEY" ]]; then
        delay=0.10
    fi

    DOWNLOADED_ACCS=()
    DOWNLOADED_PATHS=()
    local total=${#ACCESSIONS[@]}
    local i=0

    for acc in "${ACCESSIONS[@]}"; do
        i=$((i + 1))
        local out="$fasta_dir/${acc}.faa"

        if [[ -f "$out" ]]; then
            echo "  [$i/$total] $acc: already downloaded, skipping"
            DOWNLOADED_ACCS+=("$acc")
            DOWNLOADED_PATHS+=("$out")
            continue
        fi

        echo "  [$i/$total] $acc: downloading protein FASTA..."
        if download_species_faa "$acc" "$out"; then
            local count
            count=$(grep -c '^>' "$out" 2>/dev/null || echo 0)
            echo "  [$i/$total] $acc: saved $count proteins -> $out"
            DOWNLOADED_ACCS+=("$acc")
            DOWNLOADED_PATHS+=("$out")
        fi

        if [[ $i -lt $total ]]; then
            sleep "$delay"
        fi
    done

    echo ""
    echo "Downloaded ${#DOWNLOADED_ACCS[@]}/$total assemblies"
}

# ---------------------------------------------------------------------------
# Combine FASTAs
# ---------------------------------------------------------------------------
combine_fastas() {
    local output_file="$1"
    shift
    local paths=("$@")

    cat "${paths[@]}" > "$output_file"

    local count
    count=$(grep -c '^>' "$output_file" 2>/dev/null || echo 0)
    echo "Combined $count sequences into $output_file"
}

# ---------------------------------------------------------------------------
# Download ProstT5 weights
# ---------------------------------------------------------------------------
download_prostt5_weights() {
    local weights_dir="$1"
    mkdir -p "$weights_dir"
    local weights_prefix="$weights_dir/prostt5"
    local tmp_dir="$weights_dir/tmp_download"

    if [[ -e "$weights_prefix" ]]; then
        echo "ProstT5 weights already present at $weights_prefix"
        PROSTT5_WEIGHTS="$weights_prefix"
        return
    fi

    run_cmd "Downloading ProstT5 weights" \
        "$FOLDSEEK_BIN" databases ProstT5 "$weights_prefix" "$tmp_dir"
    PROSTT5_WEIGHTS="$weights_prefix"
}

# ---------------------------------------------------------------------------
# Generate 3Di codes via foldseek
# ---------------------------------------------------------------------------
generate_3di_with_foldseek() {
    local combined_fasta="$1"
    local db_dir="$2"
    mkdir -p "$db_dir"

    DB_PREFIX="$db_dir/sequenceDB"
    THREE_DI_FASTA="$db_dir/all_sequences_3di.fasta"

    # Create foldseek database with ProstT5 3Di prediction
    local cmd=("$FOLDSEEK_BIN" createdb "$combined_fasta" "$DB_PREFIX"
               --threads "$THREADS")
    if [[ -n "$PROSTT5_WEIGHTS" ]]; then
        cmd+=(--prostt5-model "$PROSTT5_WEIGHTS")
    fi
    if [[ "$USE_GPU" == true ]]; then
        cmd+=(--gpu 1)
    fi
    if [[ -n "$PROSTT5_SPLIT_LENGTH" ]]; then
        cmd+=(-prostt5SplitLength "$PROSTT5_SPLIT_LENGTH")
    fi

    run_cmd "Creating foldseek database with 3Di prediction" "${cmd[@]}"

    # Link the sequence headers to the _ss (3Di) database
    local header_db="${DB_PREFIX}_h"
    local ss_db="${DB_PREFIX}_ss"
    local ss_header="${ss_db}_h"

    run_cmd "Linking headers to 3Di database" \
        "$FOLDSEEK_BIN" lndb "$header_db" "$ss_header"

    # Extract 3Di sequences to FASTA format
    run_cmd "Extracting 3Di codes to FASTA" \
        "$FOLDSEEK_BIN" convert2fasta "$ss_db" "$THREE_DI_FASTA"

    # Count extracted sequences
    local n
    n=$(grep -c '^>' "$THREE_DI_FASTA" 2>/dev/null || echo 0)
    echo "Saved $n 3Di sequences to $THREE_DI_FASTA"
}

# ---------------------------------------------------------------------------
# Foldseek all-vs-all
# ---------------------------------------------------------------------------
run_foldseek_all_vs_all() {
    local output_dir="$1"
    mkdir -p "$output_dir"

    local result_prefix="$output_dir/result"
    RESULT_FILE="$output_dir/all_vs_all_results.tsv"
    local tmp_dir="$output_dir/tmp_search"

    # All-vs-all search
    local cmd=("$FOLDSEEK_BIN" search
               "$DB_PREFIX" "$DB_PREFIX"
               "$result_prefix" "$tmp_dir"
               --threads "$THREADS")
    if [[ ${#FOLDSEEK_ARGS[@]} -gt 0 ]]; then
        cmd+=("${FOLDSEEK_ARGS[@]}")
    fi

    run_cmd "Running foldseek all-vs-all search" "${cmd[@]}"

    # Convert to human-readable TSV
    run_cmd "Converting results to TSV" \
        "$FOLDSEEK_BIN" convertalis \
        "$DB_PREFIX" "$DB_PREFIX" \
        "$result_prefix" "$RESULT_FILE"

    echo "Results saved to $RESULT_FILE"
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
main() {
    local fasta_dir="$OUTPUT_DIR/fastas"
    local combined_fasta="$OUTPUT_DIR/combined_aa.fasta"
    local db_dir="$OUTPUT_DIR/foldseek_db"
    local metadata_file="$OUTPUT_DIR/metadata.json"
    mkdir -p "$OUTPUT_DIR"

    # Resolve foldseek binary
    resolve_foldseek_bin

    # 1. Read accessions
    read_accessions "$ACCESSION_FILE"

    # 2. Download FASTAs (or use pre-downloaded ones)
    if [[ -n "$FASTA_DIR" ]]; then
        echo ""
        echo "--- Using pre-downloaded FASTAs from $FASTA_DIR ---"
        fasta_dir="$FASTA_DIR"
        DOWNLOADED_ACCS=()
        DOWNLOADED_PATHS=()
        local total=${#ACCESSIONS[@]}
        for acc in "${ACCESSIONS[@]}"; do
            local faa_path="$fasta_dir/${acc}.faa"
            if [[ -f "$faa_path" ]]; then
                DOWNLOADED_ACCS+=("$acc")
                DOWNLOADED_PATHS+=("$faa_path")
            else
                echo "  [WARN] $acc: no .faa file found at $faa_path"
            fi
        done
        echo "Found ${#DOWNLOADED_ACCS[@]}/$total pre-downloaded FASTAs"
    else
        echo ""
        echo "--- Downloading FASTAs ---"
        batch_download "$fasta_dir"
    fi

    if [[ ${#DOWNLOADED_ACCS[@]} -eq 0 ]]; then
        echo "No sequences found. Exiting."
        exit 1
    fi

    # 3. Combine into a single FASTA
    echo ""
    echo "--- Combining sequences ---"
    combine_fastas "$combined_fasta" "${DOWNLOADED_PATHS[@]}"

    # 4. Download ProstT5 weights if not provided
    if [[ -z "$PROSTT5_WEIGHTS" ]]; then
        echo ""
        echo "--- Downloading ProstT5 weights ---"
        download_prostt5_weights "$OUTPUT_DIR/prostt5_weights"
    fi

    # 5. Generate 3Di codes via foldseek createdb
    echo ""
    echo "--- Generating 3Di codes via foldseek ---"
    generate_3di_with_foldseek "$combined_fasta" "$db_dir"

    # 6. Foldseek all-vs-all
    RESULT_FILE=""
    if [[ "$SKIP_FOLDSEEK" == true ]]; then
        echo ""
        echo "Skipping foldseek search (--skip-foldseek set)."
    else
        echo ""
        echo "--- Running foldseek all-vs-all ---"
        if ! run_foldseek_all_vs_all "$OUTPUT_DIR/foldseek_results"; then
            echo ""
            echo "Foldseek search failed. 3Di sequences are available at $THREE_DI_FASTA"
        fi
    fi

    # 7. Save metadata
    cat > "$metadata_file" <<METAEOF
{
  "accession_file": "$ACCESSION_FILE",
  "assemblies_requested": [$(printf '"%s",' "${ACCESSIONS[@]}" | sed 's/,$//')],
  "assemblies_downloaded": [$(printf '"%s",' "${DOWNLOADED_ACCS[@]}" | sed 's/,$//')],
  "num_sequences": $(grep -c '^>' "$combined_fasta" 2>/dev/null || echo 0),
  "three_di_fasta": "$THREE_DI_FASTA",
  "foldseek_results": ${RESULT_FILE:+\"$RESULT_FILE\"}${RESULT_FILE:-null}
}
METAEOF

    echo ""
    echo "Done. All outputs in $OUTPUT_DIR/"
}

main
