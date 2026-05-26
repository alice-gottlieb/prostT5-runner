#!/usr/bin/env bash
# Search extracted Pfam-match 3Di sequences against a Foldseek 3Di database.
#
# Usage:
#   ./foldseek_search_extracted_3di.sh <with_3di.domtblout|dir> <foldseek_db_dir> <output_dir> [threads] [foldseek_bin]
#
# Example:
#   ./foldseek_search_extracted_3di.sh \
#       pfam_hits2-1by1 \
#       results_fs_only_full_fs_compare/foldseek_db \
#       extracted_3di_foldseek_hits \
#       8

set -u

if [[ $# -lt 3 || $# -gt 5 ]]; then
    echo "Usage: $0 <with_3di.domtblout|dir> <foldseek_db_dir> <output_dir> [threads] [foldseek_bin]" >&2
    exit 1
fi

INPUT_PATH=$1
FOLDSEEK_DB_DIR=$2
OUTPUT_DIR=$3
THREADS=${4:-4}
FOLDSEEK_BIN=${5:-foldseek}

# Validate the few assumptions that would make the later Foldseek call cryptic.
if ! [[ "$THREADS" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: threads must be a positive integer: $THREADS" >&2
    exit 1
fi

if [[ ! -e "$INPUT_PATH" ]]; then
    echo "ERROR: input path not found: $INPUT_PATH" >&2
    exit 1
fi

if [[ ! -d "$FOLDSEEK_DB_DIR" ]]; then
    echo "ERROR: Foldseek DB directory not found: $FOLDSEEK_DB_DIR" >&2
    exit 1
fi

if ! command -v "$FOLDSEEK_BIN" >/dev/null 2>&1; then
    echo "ERROR: Foldseek binary not found on PATH: $FOLDSEEK_BIN" >&2
    exit 1
fi

TARGET_DB=""
# Infer the main Foldseek DB prefix from any DB files present. Some merged
# directories have <prefix>.dbtype, while others only have companion dbtypes
# like <prefix>_h.dbtype and <prefix>_ss_h.dbtype.
CANDIDATE_PREFIXES=()
for dbtype_file in "$FOLDSEEK_DB_DIR"/*.dbtype; do
    [[ -f "$dbtype_file" ]] || continue

    db_prefix=${dbtype_file%.dbtype}
    db_name=$(basename "$db_prefix")
    db_name=${db_name%_ss_h}
    db_name=${db_name%_ss}
    db_name=${db_name%_h}
    CANDIDATE_PREFIXES+=("$FOLDSEEK_DB_DIR/$db_name")
done

# Also inspect .index files so a missing main <prefix>.dbtype does not prevent
# detection when the usable Foldseek DB files are otherwise present.
for index_file in "$FOLDSEEK_DB_DIR"/*.index; do
    [[ -f "$index_file" ]] || continue

    db_prefix=${index_file%.index}
    db_name=$(basename "$db_prefix")
    db_name=${db_name%_ss_h}
    db_name=${db_name%_ss}
    db_name=${db_name%_h}
    CANDIDATE_PREFIXES+=("$FOLDSEEK_DB_DIR/$db_name")
done

for db_prefix in "${CANDIDATE_PREFIXES[@]}"; do
    if [[ -f "$db_prefix.index" \
        && -f "${db_prefix}_ss.index" \
        && -f "${db_prefix}_h.index" \
        && -f "${db_prefix}_ss_h.index" ]]; then
        TARGET_DB="$db_prefix"
        break
    fi
done

if [[ -z "$TARGET_DB" ]]; then
    echo "ERROR: could not infer Foldseek DB prefix from DB files in $FOLDSEEK_DB_DIR" >&2
    echo "ERROR: expected files like <prefix>.index, <prefix>_ss.index, <prefix>_h.index, and <prefix>_ss_h.index" >&2
    echo "ERROR: prefix candidates are inferred from any *.dbtype or *.index files present" >&2
    echo "Found DB-ish files:" >&2
    find "$FOLDSEEK_DB_DIR" -maxdepth 1 \( -name "*.dbtype" -o -name "*.index" \) -print | sort >&2
    exit 1
fi

echo "Detected Foldseek DB prefix: $(basename "$TARGET_DB")"

mkdir -p "$OUTPUT_DIR"

QUERY_FASTA="$OUTPUT_DIR/extracted_3di_queries.fasta"
METADATA_TSV="$OUTPUT_DIR/extracted_3di_queries.metadata.tsv"
QUERY_DB="$OUTPUT_DIR/query_db/queryDB"
RESULT_PREFIX="$OUTPUT_DIR/result"
RESULT_TSV="$OUTPUT_DIR/results.tsv"
RESULT_FULL_TSV="$OUTPUT_DIR/results.full.tsv"
RESULT_FULL_RAW_TSV="$OUTPUT_DIR/results.full.raw.tsv"
TARGET_METADATA_TSV="$OUTPUT_DIR/target_metadata.tsv"
TMP_DIR="$OUTPUT_DIR/tmp_search"

rm -rf "$OUTPUT_DIR/query_db" "$TMP_DIR"
mkdir -p "$OUTPUT_DIR/query_db"

# Accept either one augmented domtblout file or a directory of them.
if [[ -d "$INPUT_PATH" ]]; then
    mapfile -t DOMTBLOUT_FILES < <(find "$INPUT_PATH" -maxdepth 1 -type f -name "*.with_3di.domtblout" | sort)
else
    DOMTBLOUT_FILES=("$INPUT_PATH")
fi

if [[ ${#DOMTBLOUT_FILES[@]} -eq 0 ]]; then
    echo "ERROR: no *.with_3di.domtblout files found in $INPUT_PATH" >&2
    exit 1
fi

echo -e "query_id\tdomtblout_file\trow_number\ttarget_id\tpfam_name\tpfam_accession\tdomain_index\tdomain_total\tali_from\tali_to\tthree_di_sequence" > "$METADATA_TSV"
: > "$QUERY_FASTA"
: > "$QUERY_DB"
: > "$QUERY_DB.index"
: > "${QUERY_DB}_ss"
: > "${QUERY_DB}_ss.index"
: > "${QUERY_DB}_h"
: > "${QUERY_DB}_h.index"
: > "${QUERY_DB}_ss_h"
: > "${QUERY_DB}_ss_h.index"
: > "$OUTPUT_DIR/query_db/queryDB.lookup"

query_count=0
seq_offset=0
header_offset=0

# Convert each non-empty extracted 3Di sequence into one Foldseek query record.
for domtblout in "${DOMTBLOUT_FILES[@]}"; do
    echo "Reading extracted 3Di sequences from $domtblout"

    while IFS=$'\t' read -r row_number target_id pfam_name pfam_accession domain_index domain_total ali_from ali_to three_di_sequence; do
        if [[ -z "$three_di_sequence" || "$three_di_sequence" == "NA" ]]; then
            continue
        fi

        query_id="${target_id}|${pfam_accession}|domain${domain_index}of${domain_total}|ali${ali_from}-${ali_to}|row${row_number}"
        safe_query_id=${query_id//[[:space:]]/_}
        header="${safe_query_id}"$'\n'
        header_len=$((${#header} + 1))
        seq_len=${#three_di_sequence}

        # Foldseek DB files store raw sequences plus byte offsets in .index files.
        printf ">%s\n%s\n" "$safe_query_id" "$three_di_sequence" >> "$QUERY_FASTA"
        printf "%s" "$three_di_sequence" >> "$QUERY_DB"
        printf "%s" "$three_di_sequence" >> "${QUERY_DB}_ss"
        printf "%s\n\0" "$safe_query_id" >> "${QUERY_DB}_h"
        printf "%s\n\0" "$safe_query_id" >> "${QUERY_DB}_ss_h"
        printf "%d\t%d\t%d\n" "$query_count" "$seq_offset" "$seq_len" >> "$QUERY_DB.index"
        printf "%d\t%d\t%d\n" "$query_count" "$seq_offset" "$seq_len" >> "${QUERY_DB}_ss.index"
        printf "%d\t%d\t%d\n" "$query_count" "$header_offset" "$header_len" >> "${QUERY_DB}_h.index"
        printf "%d\t%d\t%d\n" "$query_count" "$header_offset" "$header_len" >> "${QUERY_DB}_ss_h.index"
        printf "%d\t%s\t0\n" "$query_count" "$safe_query_id" >> "$OUTPUT_DIR/query_db/queryDB.lookup"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$safe_query_id" "$domtblout" "$row_number" "$target_id" "$pfam_name" \
            "$pfam_accession" "$domain_index" "$domain_total" "$ali_from" "$ali_to" \
            "$three_di_sequence" >> "$METADATA_TSV"

        query_count=$((query_count + 1))
        seq_offset=$((seq_offset + seq_len))
        header_offset=$((header_offset + header_len))
    done < <(
        # The 3Di sequence is appended after a tab; the original domtblout row is
        # still space-aligned, so split that first field to recover HMMER columns.
        awk -F '\t' '
            !/^#/ && NF >= 2 && length($NF) > 0 {
                split($1, f, /[[:space:]]+/)
                print NR "\t" f[1] "\t" f[4] "\t" f[5] "\t" f[10] "\t" f[11] "\t" f[18] "\t" f[19] "\t" $NF
            }
        ' "$domtblout"
    )
done

if [[ "$query_count" -eq 0 ]]; then
    echo "ERROR: no non-empty extracted 3Di sequences found in input." >&2
    exit 1
fi

# Foldseek DB type files: 0 for sequence DB, 12 for header DB, matching DBs
# generated by Foldseek createdb for this repo's 3Di sub-databases.
printf '\000\000\000\000' > "$QUERY_DB.dbtype"
printf '\000\000\000\000' > "${QUERY_DB}_ss.dbtype"
printf '\014\000\000\000' > "${QUERY_DB}_h.dbtype"
printf '\014\000\000\000' > "${QUERY_DB}_ss_h.dbtype"

FASTAS_DIR="$(dirname "$FOLDSEEK_DB_DIR")/fastas"
echo -e "target_id\ttarget_genome\ttarget_species\ttarget_genome_species" > "$TARGET_METADATA_TSV"
# Recover target genome/species labels from the original FAA files when present.
if [[ -d "$FASTAS_DIR" ]]; then
    for faa in "$FASTAS_DIR"/*.faa; do
        [[ -f "$faa" ]] || continue
        genome=$(basename "$faa" .faa)
        awk -v genome="$genome" '
            /^>/ {
                target_id = substr($1, 2)
                species = "NA"
                if (match($0, /\[[^]]+\]/)) {
                    species = substr($0, RSTART + 1, RLENGTH - 2)
                }
                print target_id "\t" genome "\t" species "\t" genome "|" species
            }
        ' "$faa" >> "$TARGET_METADATA_TSV"
    done
    echo "Target metadata: $TARGET_METADATA_TSV"
else
    echo "WARN: target FASTA directory not found: $FASTAS_DIR"
    echo "WARN: results.full.tsv target_genome_species values will be NA"
fi

echo "Built $query_count extracted 3Di query records:"
echo "  FASTA:    $QUERY_FASTA"
echo "  Metadata: $METADATA_TSV"
echo "  Query DB: $QUERY_DB"
echo "Searching against: $TARGET_DB"

# Search all extracted 3Di domain queries against the detected Foldseek DB.
"$FOLDSEEK_BIN" search \
    "$QUERY_DB" \
    "$TARGET_DB" \
    "$RESULT_PREFIX" \
    "$TMP_DIR" \
    --threads "$THREADS"

"$FOLDSEEK_BIN" convertalis \
    "$QUERY_DB" \
    "$TARGET_DB" \
    "$RESULT_PREFIX" \
    "$RESULT_TSV" \
    --threads "$THREADS"

# Keep the raw full-format file, then add a final target_genome_species column.
"$FOLDSEEK_BIN" convertalis \
    "$QUERY_DB" \
    "$TARGET_DB" \
    "$RESULT_PREFIX" \
    "$RESULT_FULL_RAW_TSV" \
    --threads "$THREADS" \
    --format-output "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"

awk -F '\t' -v OFS='\t' '
    NR == FNR {
        if (FNR > 1) {
            target_genome_species[$1] = $4
        }
        next
    }
    {
        genome_species = ($2 in target_genome_species) ? target_genome_species[$2] : "NA"
        print $0, genome_species
    }
' "$TARGET_METADATA_TSV" "$RESULT_FULL_RAW_TSV" > "$RESULT_FULL_TSV"

echo "Done."
echo "  Default results: $RESULT_TSV"
echo "  Full results:    $RESULT_FULL_TSV"
echo "  Query metadata:  $METADATA_TSV"
echo "  Target metadata: $TARGET_METADATA_TSV"
