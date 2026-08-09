#!/bin/bash
# Create sample_name.txt and arrange files into per-sample sub-directories.
#
# Supports two input layouts inside <sample_dir>:
#   1) flat layout:    <dir>/<sample>.mt.no.softclip.bam, <dir>/<sample>.bai, ...
#                      (files are moved into a per-sample sub-directory)
#   2) nested layout:  <dir>/<sample>/<sample>.mt.no.softclip.bam, ... (already organized)
#
# Usage: bash organize_samples.sh <sample_dir>

set -euo pipefail

SAMPLE_DIR="${1:?Usage: organize_samples.sh <sample_dir>}"
cd "$SAMPLE_DIR"

if [ ! -f sample_name.txt ]; then
    # prefer the no-softclip bam names used by the feature extractor
    bams=$(find . -name "*.mt.no.softclip.bam" | sed 's#^\./##')
    if [ -z "$bams" ]; then
        bams=$(find . -name "*.mt.bam" | sed 's#^\./##')
    fi
    if [ -z "$bams" ]; then
        echo "ERROR: no *.mt.no.softclip.bam / *.mt.bam files found under $SAMPLE_DIR" >&2
        exit 1
    fi

    printf '%s\n' "$bams" \
        | sed 's#.*/##' \
        | sed 's/\.mt\.no\.softclip\.bam$//; s/\.mt\.bam$//' \
        | sort -u > sample_name.txt
    echo "Created sample_name.txt with $(wc -l < sample_name.txt) sample(s)"
else
    echo "Reusing existing sample_name.txt ($(wc -l < sample_name.txt) sample(s))"
fi

# move top-level sample files into per-sample sub-directories (nested layout is left untouched)
while IFS= read -r sample; do
    [ -z "$sample" ] && continue
    mkdir -p "$sample"
    for f in ${sample}.*; do
        if [ -f "$f" ]; then
            mv "$f" "$sample/"
        fi
    done
done < sample_name.txt

echo "organize_samples.sh finished"
