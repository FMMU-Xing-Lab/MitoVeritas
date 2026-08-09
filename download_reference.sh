#!/bin/bash
# Download the human mitochondrial reference genome (rCRS, NC_012920.1) and
# save it as reference/human_mtDNA.fasta with the sequence named "chrM".
# Usage: bash download_reference.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="$SCRIPT_DIR/reference/human_mtDNA.fasta"
TMP="$(mktemp)"

echo "Downloading NC_012920.1 (rCRS) from NCBI..."
curl -fsSL "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_012920.1&rettype=fasta&retmode=text" -o "$TMP"

# rename the header to chrM
{
    echo ">chrM"
    grep -v "^>" "$TMP"
} > "$OUT"
rm -f "$TMP"

bases=$(grep -v "^>" "$OUT" | tr -d "\n" | wc -c)
echo "Reference written to $OUT ($bases bp)"

if command -v samtools >/dev/null 2>&1; then
    samtools faidx "$OUT"
    echo "Index created: $OUT.fai"
fi
