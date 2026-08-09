#!/bin/bash
# mtDNA-ML-Predictor master pipeline
#
# End-to-end workflow:
#   mtDNApipe outputs (per-sample bam + hetro/homo mutation files)
#     -> feature extraction (R + Python)
#     -> delete 6 manual-drop features
#     -> add 9 annotation features
#     -> train XGBoost on the HCC training set and predict new mutations
#
# Usage:
#   bash run_pipeline.sh -i <sample_dir> [options]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ------------------------- defaults -------------------------
CUTOFF=2                # cutoff (%) for low-frequency variants (filter disabled by default;
                        # set APPLY_CUTOFF=1 to enable it)
VAF=0.1                 # vaf tag used in hetro/homo file names
THREADS=4               # parallel feature-extraction jobs
REFERENCE="${MTDNAPIPE_REFERENCE:-$SCRIPT_DIR/reference/human_mtDNA.fasta}"
TRAIN="${MTDNAPIPE_TRAIN:-$SCRIPT_DIR/data/HCC_training_all.txt}"
OUTPUT_DIR=""
SAMPLE_NAME="prediction"
SKIP_CHECK=0

usage() {
    cat <<EOF
Usage: bash run_pipeline.sh -i <sample_dir> [options]

Required:
  -i DIR   analysis folder containing per-sample sub-folders (or flat sample files)

Options:
  -o DIR   output directory for final results (default: <sample_dir>/outputs)
  -c INT   cutoff (%) for low-frequency variants (default: $CUTOFF; the filter is
           disabled by default, enable it with the environment variable APPLY_CUTOFF=1)
  -v FLOAT vaf tag used in the hetro/homo file names (default: $VAF)
  -t INT   number of parallel feature-extraction jobs (default: $THREADS)
  -r FILE  mitochondrial reference fasta (default: $REFERENCE)
  -T FILE  training set for prediction (default: $TRAIN)
  -s NAME  sample/cohort name used in the prediction log (default: $SAMPLE_NAME)
  -k       skip the environment check
  -h       show this help

Environment variables:
  MTDNAPIPE_REFERENCE   alternative reference fasta path
  MTDNAPIPE_TRAIN       alternative training set path
  APPLY_CUTOFF=1        apply the cutoff filter when preparing candidate variants
EOF
}

while getopts "i:o:c:v:t:r:T:s:kh" opt; do
    case "$opt" in
        i) SAMPLE_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        c) CUTOFF="$OPTARG" ;;
        v) VAF="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        r) REFERENCE="$OPTARG" ;;
        T) TRAIN="$OPTARG" ;;
        s) SAMPLE_NAME="$OPTARG" ;;
        k) SKIP_CHECK=1 ;;
        h) usage; exit 0 ;;
        *) usage; exit 1 ;;
    esac
done

if [ -z "${SAMPLE_DIR:-}" ]; then
    usage
    exit 1
fi

INPUT_DIR="$(cd "$SAMPLE_DIR" && pwd)"
OUTPUT_DIR="${OUTPUT_DIR:-$INPUT_DIR/outputs}"
mkdir -p "$OUTPUT_DIR"

echo "=================================================="
echo " mtDNA-ML-Predictor"
echo "=================================================="
echo "Input dir    : $INPUT_DIR"
echo "Output dir   : $OUTPUT_DIR"
echo "Cutoff       : $CUTOFF"
echo "VAF tag      : $VAF"
echo "Threads      : $THREADS"
echo "Reference    : $REFERENCE"
echo "Training set : $TRAIN"
echo "Sample name  : $SAMPLE_NAME"

# ------------------------- 0. environment check -------------------------
if [ "$SKIP_CHECK" -eq 0 ]; then
    echo
    echo ">>> Step 0: checking environment"
    bash "$SCRIPT_DIR/check_environment.sh"
fi

[ -f "$REFERENCE" ] || { echo "ERROR: reference fasta not found: $REFERENCE" >&2; exit 1; }
[ -f "$TRAIN" ] || { echo "ERROR: training set not found: $TRAIN" >&2; exit 1; }

export MTDNAPIPE_REFERENCE="$REFERENCE"

if [ ! -f "$REFERENCE.fai" ]; then
    echo ">>> Creating reference index ($REFERENCE.fai)"
    samtools faidx "$REFERENCE"
fi

# ------------------------- 1. feature extraction -------------------------
echo
echo ">>> Step 1/4: feature extraction from mtDNApipe outputs"
bash "$SCRIPT_DIR/scripts/get_feature_pipeline.sh" \
    "$INPUT_DIR" "$CUTOFF" "$THREADS" "$VAF" "$REFERENCE"

FEATURE_RAW="$INPUT_DIR/Total_output_feature_${VAF}.txt"
[ -f "$FEATURE_RAW" ] || { echo "ERROR: aggregated feature file not found: $FEATURE_RAW" >&2; exit 1; }

# ------------------------- 2. delete features -------------------------
echo
echo ">>> Step 2/4: deleting the 6 manual-drop features"
python3 "$SCRIPT_DIR/scripts/delete_features.py" \
    --input "$FEATURE_RAW" \
    --output "$OUTPUT_DIR/feature.txt"

# ------------------------- 3. add annotation features -------------------------
echo
echo ">>> Step 3/4: adding annotation features"
python3 "$SCRIPT_DIR/scripts/add_features.py" \
    --input "$OUTPUT_DIR/feature.txt" \
    --output "$OUTPUT_DIR/feature_add10.txt" \
    --data-dir "$SCRIPT_DIR/data"

# ------------------------- 4. predict -------------------------
echo
echo ">>> Step 4/4: training and predicting"
python3 "$SCRIPT_DIR/scripts/predict.py" \
    --train "$TRAIN" \
    --predict "$OUTPUT_DIR/feature_add10.txt" \
    --output "$OUTPUT_DIR/${SAMPLE_NAME}-pred-results.txt" \
    --sample-name "$SAMPLE_NAME"

echo
echo "=================================================="
echo "Pipeline finished successfully!"
echo "  feature.txt           : $OUTPUT_DIR/feature.txt"
echo "  feature_add10.txt     : $OUTPUT_DIR/feature_add10.txt"
echo "  prediction results    : $OUTPUT_DIR/${SAMPLE_NAME}-pred-results.txt"
echo "=================================================="
