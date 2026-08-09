#!/bin/bash
# Extract ML features from mtDNApipe outputs (portable version of get_mtDNApipe_feature.sh).
#
# Usage: bash get_feature_pipeline.sh <sample_dir> <cutoff> <threads> <vaf> [reference.fasta]
#
#   sample_dir : analysis folder containing sample_name.txt (or flat sample files)
#   cutoff     : cutoff (%) for low-frequency variants (filter disabled by default,
#                set APPLY_CUTOFF=1 to enable)
#   threads    : number of parallel feature-extraction jobs
#   vaf        : mutation-frequency threshold tag used in hetro/homo file names (e.g. 0.1)
#   reference  : mitochondrial reference fasta (default: $MTDNAPIPE_REFERENCE)

set -euo pipefail

SAMPLE_DIR="$1"
CUTOFF="${2:-2}"
THREADS="${3:-4}"
VAF="${4:-0.1}"
REFERENCE="${5:-${MTDNAPIPE_REFERENCE:-}}"

if [ -z "$REFERENCE" ]; then
    echo "ERROR: reference fasta not provided (argument 5 or MTDNAPIPE_REFERENCE)" >&2
    exit 1
fi
export MTDNAPIPE_REFERENCE="$REFERENCE"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

[ -d "$SAMPLE_DIR" ] || { echo "ERROR: sample dir not found: $SAMPLE_DIR" >&2; exit 1; }
[ -f "$REFERENCE" ] || { echo "ERROR: reference fasta not found: $REFERENCE" >&2; exit 1; }
if [ ! -f "$REFERENCE.fai" ]; then
    echo "Creating reference index: $REFERENCE.fai"
    samtools faidx "$REFERENCE"
fi

start=$(date +%s)
echo "Start time   : $(date)"
echo "Sample dir   : $SAMPLE_DIR"
echo "Cutoff       : $CUTOFF"
echo "Threads      : $THREADS"
echo "VAF tag      : $VAF"
echo "Reference    : $REFERENCE"
echo "Scripts dir  : $SCRIPT_DIR"

cd "$SAMPLE_DIR"

# remove previous intermediate outputs for this vaf
find . -maxdepth 2 -name "output_true_${VAF}" -exec rm -f {} \;
find . -maxdepth 2 -name "output_feature_${VAF}*" -exec rm -f {} \;

# ---- (1) organize samples ----
bash "$SCRIPT_DIR/organize_samples.sh" "$SAMPLE_DIR"

# ---- (2) prepare candidate variant sites (output_true_<vaf>) ----
echo "get_output_true_file.R starts at $(date)"
Rscript "$SCRIPT_DIR/get_output_true_file.R" "$SAMPLE_DIR" "$CUTOFF" "$VAF"

# ---- (3) extract features per sample in parallel ----
if [ ! -s sample_name.txt ]; then
    echo "ERROR: sample_name.txt is empty" >&2
    exit 1
fi
echo "get_feature_mtDNApipe.py starts at $(date) (threads=$THREADS)"
cat sample_name.txt | xargs -P "$THREADS" -I{} python3 "$SCRIPT_DIR/get_feature_mtDNApipe.py" "$SAMPLE_DIR" "{}" "$REFERENCE"

# ---- (4) aggregate all per-sample features ----
echo "get_total_output_features.R starts at $(date)"
Rscript "$SCRIPT_DIR/get_total_output_features.R" "$SAMPLE_DIR" "$CUTOFF"

end=$(date +%s)
elapsed=$((end - start))
echo "Feature extraction finished at $(date), elapsed $((elapsed / 60)) min $((elapsed % 60)) s"
echo "Aggregated feature file: $SAMPLE_DIR/Total_output_feature_${VAF}.txt"
