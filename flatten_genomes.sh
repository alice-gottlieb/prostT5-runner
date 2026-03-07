#!/bin/bash
# Flatten downloaded NCBI genomes: rename protein.faa to [accession].faa
# and move to a single top-level directory, then clean up subdirectories.
#
# Usage: bash flatten_genomes.sh /path/to/ncbi_genomes

GENOMES_DIR="${1:?Usage: bash flatten_genomes.sh $SCRATCH/ncbi_genomes}"

# Find all protein.faa files and move them with accession-based names
find "$GENOMES_DIR" -path "*/GCF_*/protein.faa" -o -path "*/GCA_*/protein.faa" | while read -r faa; do
    accession=$(basename "$(dirname "$faa")")
    cp "$faa" "${GENOMES_DIR}/${accession}.faa"
    echo "Copied ${accession}.faa"
done

# Remove all task subdirectories
# rm -rf "${GENOMES_DIR}"/task_*

echo "Done. Flattened files in ${GENOMES_DIR}/"
