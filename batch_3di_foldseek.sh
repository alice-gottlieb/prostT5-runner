#!/usr/bin/env bash
# Batch download whole-species protein FASTA (.faa) files from NCBI Datasets,
# convert all sequences to 3Di codes using foldseek's integrated ProstT5, then
# run an all-vs-all foldseek search across all sequences.
#
# Input:  a text file with one NCBI genome assembly accession (GCF_* or GCA_*)
#         per line
# Output: per-assembly .faa files, a combined 3Di FASTA, and foldseek results
#
# Requires: foldseek (v9+ with ProstT5 support), curl, unzip, jq
#
# Example usages:
#
#   # Basic usage (auto-downloads ProstT5 weights, uses CPU)
#   bash batch_3di_foldseek.sh assemblies.txt
#
#   # Custom output directory
#   bash batch_3di_foldseek.sh assemblies.txt -o my_results
#
#   # With GPU acceleration
#   bash batch_3di_foldseek.sh assemblies.txt --gpu
#
#   # Multi-threaded with GPU
#   bash batch_3di_foldseek.sh assemblies.txt --gpu --threads 8
#
#   # With pre-downloaded ProstT5 weights
#   bash batch_3di_foldseek.sh assemblies.txt --prostt5-weights /path/to/prostt5
#
#   # With a specific foldseek binary
#   bash batch_3di_foldseek.sh assemblies.txt --foldseek-path /opt/foldseek/bin/foldseek
#
#   # With NCBI API key (faster downloads, 10 req/s vs 3 req/s)
#   bash batch_3di_foldseek.sh assemblies.txt --ncbi-api-key YOUR_KEY
#
#   # Only generate 3Di codes, skip the all-vs-all search
#   bash batch_3di_foldseek.sh assemblies.txt --skip-foldseek
#
#   # Pass extra arguments to foldseek search (e.g. e-value cutoff)
#   bash batch_3di_foldseek.sh assemblies.txt --foldseek-args "-e 0.001"
#
#   # Full example with multiple options
#   bash batch_3di_foldseek.sh assemblies.txt \
#       -o results --gpu --threads 8 \
#       --ncbi-api-key YOUR_KEY \
#       --foldseek-args "-e 0.001 --max-seqs 1000"

set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

OUTPUT_DIR="batch_output"
FASTA_DIR=""
NCBI_API_KEY=""
FOLDSEEK_PATH=""
PROSTT5_WEIGHTS=""
USE_GPU=0
THREADS=1
SKIP_FOLDSEEK=0
FOLDSEEK_ARGS=""
PROSTT5_SPLIT_LENGTH=""

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

usage() {
    sed -n '/^# Example usages/,/^[^#]/{ /^#/p }' "$0" | sed 's/^# \{0,2\}//'
    echo ""
    echo "Usage: $0 <accession_file> [options]"
    echo ""
    echo "Options:"
    echo "  -o, --output-dir DIR        Output directory (default: batch_output)"
    echo "  --fasta-dir DIR             Use pre-downloaded .faa files, skip NCBI download"
    echo "  --ncbi-api-key KEY          NCBI API key for higher rate limits"
    echo "  --foldseek-path PATH        Path to foldseek binary or install dir"
    echo "  --prostt5-weights PATH      Path prefix for pre-downloaded ProstT5 weights"
    echo "  --gpu                       Use GPU for ProstT5 inference"
    echo "  --prostt5-split-length N    Max sequence length before splitting"
    echo "  -t, --threads N             Number of threads (default: 1)"
    echo "  --skip-foldseek             Stop after 3Di generation, skip search"
    echo "  --foldseek-args \"ARGS\"      Extra args forwarded to foldseek search"
    echo "  -h, --help                  Show this help"
    exit 0
}

if [[ $# -lt 1 ]]; then
    echo "Error: accession_file is required." >&2
    usage
fi

ACCESSION_FILE="$1"
shift

while [[ $# -gt 0 ]]; do
    case "$1" in
        -o|--output-dir)       OUTPUT_DIR="$2";           shift 2 ;;
        --fasta-dir)           FASTA_DIR="$2";             shift 2 ;;
        --ncbi-api-key)        NCBI_API_KEY="$2";          shift 2 ;;
        --foldseek-path)       FOLDSEEK_PATH="$2";         shift 2 ;;
        --prostt5-weights)     PROSTT5_WEIGHTS="$2";       shift 2 ;;
        --gpu)                 USE_GPU=1;                  shift   ;;
        --prostt5-split-length) PROSTT5_SPLIT_LENGTH="$2"; shift 2 ;;
        -t|--threads)          THREADS="$2";               shift 2 ;;
        --skip-foldseek)       SKIP_FOLDSEEK=1;            shift   ;;
        --foldseek-args)       FOLDSEEK_ARGS="$2";         shift 2 ;;
        -h|--help)             usage ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

log()  { echo "$*"; }
warn() { echo "  [WARN] $*"; }

run_cmd() {
    local description="$1"; shift
    log ""
    log "${description}: $*"
    "$@"
}

# ---------------------------------------------------------------------------
# Resolve foldseek binary
# ---------------------------------------------------------------------------

resolve_foldseek_bin() {
    if [[ -n "$FOLDSEEK_PATH" ]]; then
        if [[ -d "$FOLDSEEK_PATH" ]]; then
            local bin="${FOLDSEEK_PATH}/bin/foldseek"
            [[ -f "$bin" ]] || bin="${FOLDSEEK_PATH}/foldseek"
        else
            local bin="$FOLDSEEK_PATH"
        fi
        if [[ ! -f "$bin" ]]; then
            echo "Error: foldseek binary not found at ${bin}" >&2
            exit 1
        fi
        echo "$bin"
    else
        echo "foldseek"
    fi
}

# ---------------------------------------------------------------------------
# Read accessions
# ---------------------------------------------------------------------------

read_accessions() {
    local file="$1"
    # Filter blank lines and comments
    grep -v '^\s*$' "$file" | grep -v '^\s*#' || true
}

# ---------------------------------------------------------------------------
# Download a single .faa from NCBI Datasets
# ---------------------------------------------------------------------------

download_species_faa() {
    local accession="$1"
    local output_path="$2"

    local url="https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/${accession}/download"
    local tmp_zip
    tmp_zip=$(mktemp /tmp/ncbi_XXXXXX.zip)

    local curl_args=(-fsSL --max-time 120 -o "$tmp_zip"
                     -G "$url"
                     --data-urlencode "include_annotation_type=PROT_FASTA"
                     -H "Accept: application/zip")
    if [[ -n "$NCBI_API_KEY" ]]; then
        curl_args+=(-H "api-key: ${NCBI_API_KEY}")
    fi

    if ! curl "${curl_args[@]}"; then
        warn "${accession}: download failed"
        rm -f "$tmp_zip"
        return 1
    fi

    # Extract protein.faa from the zip
    local faa_entry
    faa_entry=$(unzip -l "$tmp_zip" 2>/dev/null | awk '/protein\.faa$/{print $NF}' | head -1)
    if [[ -z "$faa_entry" ]]; then
        warn "${accession}: no protein.faa found in downloaded zip"
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
    local -a accessions=("$@")
    local total=${#accessions[@]}
    local delay=0.34
    [[ -n "$NCBI_API_KEY" ]] && delay=0.10

    local i=0
    for acc in "${accessions[@]}"; do
        (( i++ )) || true
        local out="${FASTA_DIR_ACTUAL}/${acc}.faa"
        if [[ -f "$out" ]]; then
            log "  [${i}/${total}] ${acc}: already downloaded, skipping"
            continue
        fi
        log "  [${i}/${total}] ${acc}: downloading protein FASTA..."
        if download_species_faa "$acc" "$out"; then
            local count
            count=$(grep -c '^>' "$out" 2>/dev/null || echo 0)
            log "  [${i}/${total}] ${acc}: saved ${count} proteins -> ${out}"
        else
            log "  [${i}/${total}] ${acc}: skipped (download failed)"
        fi
        [[ $i -lt $total ]] && sleep "$delay"
    done
}

# ---------------------------------------------------------------------------
# Combine FASTAs
# ---------------------------------------------------------------------------

combine_fastas() {
    local combined="$1"; shift
    local -a faa_files=("$@")
    : > "$combined"
    for f in "${faa_files[@]}"; do
        cat "$f" >> "$combined"
    done
    local count
    count=$(grep -c '^>' "$combined" 2>/dev/null || echo 0)
    log "Combined ${count} sequences into ${combined}"
}

# ---------------------------------------------------------------------------
# Download ProstT5 weights
# ---------------------------------------------------------------------------

download_prostt5_weights() {
    local foldseek_bin="$1"
    local weights_dir="$2"
    mkdir -p "$weights_dir"
    local weights_prefix="${weights_dir}/prostt5"
    local tmp_dir="${weights_dir}/tmp_download"

    if [[ -f "${weights_prefix}" || -d "${weights_prefix}" ]]; then
        log "ProstT5 weights already present at ${weights_prefix}"
        echo "$weights_prefix"
        return
    fi

    run_cmd "Downloading ProstT5 weights" \
        "$foldseek_bin" databases ProstT5 "$weights_prefix" "$tmp_dir"
    echo "$weights_prefix"
}

# ---------------------------------------------------------------------------
# Generate 3Di via foldseek createdb
# ---------------------------------------------------------------------------

generate_3di_with_foldseek() {
    local combined_fasta="$1"
    local db_dir="$2"
    local foldseek_bin="$3"
    local weights="$4"

    mkdir -p "$db_dir"
    local db_prefix="${db_dir}/sequenceDB"
    local three_di_fasta="${db_dir}/all_sequences_3di.fasta"

    local cmd=("$foldseek_bin" createdb "$combined_fasta" "$db_prefix"
                "--threads" "$THREADS")
    [[ -n "$weights" ]]              && cmd+=("--prostt5-model" "$weights")
    [[ "$USE_GPU" -eq 1 ]]           && cmd+=("--gpu" "1")
    [[ -n "$PROSTT5_SPLIT_LENGTH" ]] && cmd+=("-prostt5SplitLength" "$PROSTT5_SPLIT_LENGTH")

    run_cmd "Creating foldseek database with 3Di prediction" "${cmd[@]}"

    local header_db="${db_prefix}_h"
    local ss_db="${db_prefix}_ss"
    local ss_header="${ss_db}_h"

    run_cmd "Linking headers to 3Di database" \
        "$foldseek_bin" lndb "$header_db" "$ss_header"

    run_cmd "Extracting 3Di codes to FASTA" \
        "$foldseek_bin" convert2fasta "$ss_db" "$three_di_fasta"

    local n
    n=$(grep -c '^>' "$three_di_fasta" 2>/dev/null || echo 0)
    log "Saved ${n} 3Di sequences to ${three_di_fasta}"

    # Export for caller
    DB_PREFIX_OUT="$db_prefix"
    THREE_DI_FASTA_OUT="$three_di_fasta"
}

# ---------------------------------------------------------------------------
# Foldseek all-vs-all search
# ---------------------------------------------------------------------------

run_foldseek_all_vs_all() {
    local db_prefix="$1"
    local output_dir="$2"
    local foldseek_bin="$3"

    mkdir -p "$output_dir"
    local result_prefix="${output_dir}/result"
    local result_file="${output_dir}/all_vs_all_results.tsv"
    local tmp_dir="${output_dir}/tmp_search"

    # shellcheck disable=SC2206
    local extra=($FOLDSEEK_ARGS)

    run_cmd "Running foldseek all-vs-all search" \
        "$foldseek_bin" search \
        "$db_prefix" "$db_prefix" \
        "$result_prefix" "$tmp_dir" \
        "--threads" "$THREADS" \
        "${extra[@]+"${extra[@]}"}"

    run_cmd "Converting results to TSV" \
        "$foldseek_bin" convertalis \
        "$db_prefix" "$db_prefix" \
        "$result_prefix" "$result_file"

    log "Results saved to ${result_file}"
    RESULT_FILE_OUT="$result_file"
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

main() {
    # Validate input file
    if [[ ! -f "$ACCESSION_FILE" ]]; then
        echo "Error: accession file not found: ${ACCESSION_FILE}" >&2
        exit 1
    fi

    local out="$OUTPUT_DIR"
    local combined_fasta="${out}/combined_aa.fasta"
    local db_dir="${out}/foldseek_db"
    local metadata_file="${out}/metadata.json"
    mkdir -p "$out"

    local foldseek_bin
    foldseek_bin=$(resolve_foldseek_bin)

    # 1. Read accessions
    mapfile -t accessions < <(read_accessions "$ACCESSION_FILE")
    local total=${#accessions[@]}
    log "Loaded ${total} assembly accessions from ${ACCESSION_FILE}"

    if [[ $total -eq 0 ]]; then
        echo "No accessions found in ${ACCESSION_FILE}. Exiting." >&2
        exit 1
    fi

    # 2. Locate or download FASTAs
    if [[ -n "$FASTA_DIR" ]]; then
        log ""
        log "--- Using pre-downloaded FASTAs from ${FASTA_DIR} ---"
        FASTA_DIR_ACTUAL="$FASTA_DIR"
    else
        FASTA_DIR_ACTUAL="${out}/fastas"
        mkdir -p "$FASTA_DIR_ACTUAL"
        log ""
        log "--- Downloading FASTAs ---"
        batch_download "${accessions[@]}"
    fi

    # Collect successfully present .faa files
    local -a faa_files=()
    local -a downloaded_accs=()
    for acc in "${accessions[@]}"; do
        local faa="${FASTA_DIR_ACTUAL}/${acc}.faa"
        if [[ -f "$faa" ]]; then
            faa_files+=("$faa")
            downloaded_accs+=("$acc")
        else
            warn "${acc}: no .faa file found, skipping"
        fi
    done

    log ""
    log "Downloaded/found ${#faa_files[@]}/${total} assemblies"

    if [[ ${#faa_files[@]} -eq 0 ]]; then
        echo "No sequences found. Exiting." >&2
        exit 1
    fi

    # 3+4. Combine FASTAs
    log ""
    log "--- Combining sequences ---"
    combine_fastas "$combined_fasta" "${faa_files[@]}"

    # 5. Resolve ProstT5 weights
    if [[ -z "$PROSTT5_WEIGHTS" ]]; then
        log ""
        log "--- Downloading ProstT5 weights ---"
        PROSTT5_WEIGHTS=$(download_prostt5_weights "$foldseek_bin" "${out}/prostt5_weights")
    fi

    # 6. Generate 3Di codes
    log ""
    log "--- Generating 3Di codes via foldseek ---"
    DB_PREFIX_OUT=""
    THREE_DI_FASTA_OUT=""
    generate_3di_with_foldseek \
        "$combined_fasta" "$db_dir" "$foldseek_bin" "$PROSTT5_WEIGHTS"

    local db_prefix="$DB_PREFIX_OUT"
    local three_di_fasta="$THREE_DI_FASTA_OUT"

    # 7. All-vs-all search
    local result_file="null"
    if [[ "$SKIP_FOLDSEEK" -eq 1 ]]; then
        log ""
        log "Skipping foldseek search (--skip-foldseek set)."
    else
        log ""
        log "--- Running foldseek all-vs-all ---"
        RESULT_FILE_OUT=""
        if run_foldseek_all_vs_all "$db_prefix" "${out}/foldseek_results" "$foldseek_bin"; then
            result_file="\"${RESULT_FILE_OUT}\""
        else
            log "Foldseek search failed. 3Di sequences are available at ${three_di_fasta}"
        fi
    fi

    # 8. Save metadata JSON
    local downloaded_json
    downloaded_json=$(printf '"%s",' "${downloaded_accs[@]}")
    downloaded_json="[${downloaded_json%,}]"

    local accessions_json
    accessions_json=$(printf '"%s",' "${accessions[@]}")
    accessions_json="[${accessions_json%,}]"

    local num_sequences
    num_sequences=$(grep -c '^>' "$combined_fasta" 2>/dev/null || echo 0)

    cat > "$metadata_file" <<EOF
{
  "accession_file": "${ACCESSION_FILE}",
  "assemblies_requested": ${accessions_json},
  "assemblies_downloaded": ${downloaded_json},
  "num_sequences": ${num_sequences},
  "three_di_fasta": "${three_di_fasta}",
  "foldseek_results": ${result_file}
}
EOF

    log ""
    log "Done. All outputs in ${out}/"
}

main
