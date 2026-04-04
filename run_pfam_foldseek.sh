#!/bin/bash
# ---------------------------------------------------------------------------
# run_pfam_foldseek.sh
#
# Convert per-family Pfam amino acid FASTAs to 3Di codes using foldseek's
# integrated ProstT5 model.
#
# Input:  Directory of per-family FASTAs produced by download_pfam.py
#           pfam_fastas/
#             PF00001/PF00001.fasta
#             PF00002/PF00002.fasta
#             ...
#
# Output: Per-family 3Di FASTA files in $SCRATCH/pfam-3dis/
#           $SCRATCH/pfam-3dis/PF00001/PF00001_3di.fasta
#           $SCRATCH/pfam-3dis/PF00002/PF00002_3di.fasta
#           ...
#
# Usage:
#   # Basic (CPU, auto-download ProstT5 weights)
#   bash run_pfam_foldseek.sh pfam_fastas
#
#   # With GPU and custom foldseek binary
#   bash run_pfam_foldseek.sh pfam_fastas --gpu --foldseek /path/to/foldseek
#
#   # With pre-downloaded ProstT5 weights
#   bash run_pfam_foldseek.sh pfam_fastas --prostt5-weights /path/to/prostt5
#
#   # Process a subset (e.g., for array jobs)
#   bash run_pfam_foldseek.sh pfam_fastas --start 0 --count 1000
#
#   # Multi-threaded
#   bash run_pfam_foldseek.sh pfam_fastas --threads 4
# ---------------------------------------------------------------------------
set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
INPUT_DIR=""
OUTPUT_BASE="${SCRATCH}/pfam-3dis"
FOLDSEEK_BIN="foldseek"
PROSTT5_WEIGHTS=""
GPU_FLAG=""
THREADS=1
START=0
COUNT=0  # 0 = process all
LOG_DIR=""

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
usage() {
    echo "Usage: $0 <input_dir> [options]"
    echo ""
    echo "Arguments:"
    echo "  input_dir              Directory with per-family FASTAs from download_pfam.py"
    echo ""
    echo "Options:"
    echo "  --output-dir DIR       Output directory (default: \$SCRATCH/pfam-3dis)"
    echo "  --foldseek PATH        Path to foldseek binary (default: foldseek)"
    echo "  --prostt5-weights PATH Path prefix for ProstT5 weights (.gguf)"
    echo "  --gpu                  Use GPU for ProstT5 inference"
    echo "  --threads N            Number of threads (default: 1)"
    echo "  --start N              Start at family index N (for chunked processing)"
    echo "  --count N              Process N families (0 = all, default: 0)"
    echo "  --log-dir DIR          Directory for per-family log files"
    echo "  -h, --help             Show this help message"
    exit 1
}

if [ $# -lt 1 ]; then
    usage
fi

INPUT_DIR="$1"
shift

while [ $# -gt 0 ]; do
    case "$1" in
        --output-dir)   OUTPUT_BASE="$2"; shift 2 ;;
        --foldseek)     FOLDSEEK_BIN="$2"; shift 2 ;;
        --prostt5-weights) PROSTT5_WEIGHTS="$2"; shift 2 ;;
        --gpu)          GPU_FLAG="--gpu 1"; shift ;;
        --threads)      THREADS="$2"; shift 2 ;;
        --start)        START="$2"; shift 2 ;;
        --count)        COUNT="$2"; shift 2 ;;
        --log-dir)      LOG_DIR="$2"; shift 2 ;;
        -h|--help)      usage ;;
        *)              echo "Unknown option: $1"; usage ;;
    esac
done

# ---------------------------------------------------------------------------
# Validate
# ---------------------------------------------------------------------------
if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: Input directory not found: $INPUT_DIR"
    exit 1
fi

if ! command -v "$FOLDSEEK_BIN" &>/dev/null && [ ! -x "$FOLDSEEK_BIN" ]; then
    echo "ERROR: foldseek binary not found: $FOLDSEEK_BIN"
    echo "  Install foldseek or provide --foldseek /path/to/foldseek"
    exit 1
fi

if [ -z "$SCRATCH" ] && [ "$OUTPUT_BASE" = "\$SCRATCH/pfam-3dis" ]; then
    echo "ERROR: \$SCRATCH is not set. Either set it or use --output-dir."
    exit 1
fi

mkdir -p "$OUTPUT_BASE"

if [ -n "$LOG_DIR" ]; then
    mkdir -p "$LOG_DIR"
fi

# Build prostt5 model flag
PROSTT5_FLAG=""
if [ -n "$PROSTT5_WEIGHTS" ]; then
    PROSTT5_FLAG="--prostt5-model $PROSTT5_WEIGHTS"
fi

# ---------------------------------------------------------------------------
# Collect family directories
# ---------------------------------------------------------------------------
echo "=== Pfam 3Di Pipeline ==="
echo "  Input:       $INPUT_DIR"
echo "  Output:      $OUTPUT_BASE"
echo "  Foldseek:    $FOLDSEEK_BIN"
echo "  ProstT5:     ${PROSTT5_WEIGHTS:-auto}"
echo "  GPU:         ${GPU_FLAG:-no}"
echo "  Threads:     $THREADS"
echo ""

# Get sorted list of family directories that contain a .fasta file
FAMILIES=()
while IFS= read -r dir; do
    fam=$(basename "$dir")
    if [ -f "$dir/${fam}.fasta" ]; then
        FAMILIES+=("$fam")
    fi
done < <(find "$INPUT_DIR" -mindepth 1 -maxdepth 1 -type d | sort)

TOTAL=${#FAMILIES[@]}
echo "Found $TOTAL Pfam families in $INPUT_DIR"

if [ "$TOTAL" -eq 0 ]; then
    echo "No families found. Exiting."
    exit 0
fi

# Apply start/count slicing
if [ "$START" -gt 0 ] || [ "$COUNT" -gt 0 ]; then
    END=$TOTAL
    if [ "$COUNT" -gt 0 ]; then
        END=$(( START + COUNT ))
        if [ "$END" -gt "$TOTAL" ]; then
            END=$TOTAL
        fi
    fi
    FAMILIES=("${FAMILIES[@]:$START:$((END - START))}")
    echo "Processing families $START to $((END - 1)) (${#FAMILIES[@]} families)"
fi

# ---------------------------------------------------------------------------
# Process each family
# ---------------------------------------------------------------------------
PROCESSED=0
SKIPPED=0
FAILED=0

for FAM in "${FAMILIES[@]}"; do
    PROCESSED=$((PROCESSED + 1))
    INPUT_FASTA="$INPUT_DIR/$FAM/${FAM}.fasta"
    FAM_OUTPUT="$OUTPUT_BASE/$FAM"
    OUTPUT_3DI="$FAM_OUTPUT/${FAM}_3di.fasta"

    # Skip if 3Di output already exists
    if [ -f "$OUTPUT_3DI" ]; then
        SKIPPED=$((SKIPPED + 1))
        if [ $((PROCESSED % 500)) -eq 0 ]; then
            echo "  [$PROCESSED/${#FAMILIES[@]}] Progress update (skipped $SKIPPED so far)..."
        fi
        continue
    fi

    echo "  [$PROCESSED/${#FAMILIES[@]}] Processing $FAM..."

    mkdir -p "$FAM_OUTPUT"

    # Temp directory for foldseek database files
    DB_DIR="$FAM_OUTPUT/foldseek_db"
    DB_PREFIX="$DB_DIR/sequenceDB"
    TMP_DIR="$DB_DIR/tmp"
    mkdir -p "$DB_DIR" "$TMP_DIR"

    # Redirect logs if log dir is set
    if [ -n "$LOG_DIR" ]; then
        LOG_FILE="$LOG_DIR/${FAM}.log"
        exec 3>"$LOG_FILE"
    else
        exec 3>/dev/null
    fi

    # Step 1: Create foldseek database with ProstT5 3Di prediction
    if ! "$FOLDSEEK_BIN" createdb \
        "$INPUT_FASTA" "$DB_PREFIX" \
        --threads "$THREADS" \
        $PROSTT5_FLAG \
        $GPU_FLAG \
        >&3 2>&1; then
        echo "    FAILED: foldseek createdb for $FAM"
        FAILED=$((FAILED + 1))
        rm -rf "$DB_DIR"
        exec 3>&-
        continue
    fi

    # Step 2: Link headers to 3Di (_ss) database
    if ! "$FOLDSEEK_BIN" lndb \
        "${DB_PREFIX}_h" "${DB_PREFIX}_ss_h" \
        >&3 2>&1; then
        echo "    FAILED: foldseek lndb for $FAM"
        FAILED=$((FAILED + 1))
        rm -rf "$DB_DIR"
        exec 3>&-
        continue
    fi

    # Step 3: Extract 3Di sequences to FASTA
    if ! "$FOLDSEEK_BIN" convert2fasta \
        "${DB_PREFIX}_ss" "$OUTPUT_3DI" \
        >&3 2>&1; then
        echo "    FAILED: foldseek convert2fasta for $FAM"
        FAILED=$((FAILED + 1))
        rm -rf "$DB_DIR"
        exec 3>&-
        continue
    fi

    exec 3>&-

    # Clean up intermediate database files to save disk space
    rm -rf "$DB_DIR"

    # Report sequence count
    N_SEQS=$(grep -c '^>' "$OUTPUT_3DI" 2>/dev/null || echo 0)
    if [ $((PROCESSED % 100)) -eq 0 ]; then
        echo "    -> $N_SEQS 3Di sequences written to $OUTPUT_3DI"
        echo "  [$PROCESSED/${#FAMILIES[@]}] Progress: $SKIPPED skipped, $FAILED failed"
    fi
done

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo ""
echo "=== Complete ==="
echo "  Total families:  ${#FAMILIES[@]}"
echo "  Processed:       $((PROCESSED - SKIPPED - FAILED))"
echo "  Skipped (exist): $SKIPPED"
echo "  Failed:          $FAILED"
echo "  Output:          $OUTPUT_BASE/[PFxxxxx]/"
echo ""

if [ "$FAILED" -gt 0 ]; then
    echo "WARNING: $FAILED families failed. Check logs for details."
    exit 1
fi
