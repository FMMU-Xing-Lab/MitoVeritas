#!/bin/bash
# Check that the current environment can run the mtDNA-ML-Predictor pipeline.
# Usage: bash check_environment.sh

set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PASS=0
FAIL=0

check_cmd() {
    if command -v "$1" >/dev/null 2>&1; then
        echo "[OK]   $1 -> $(command -v "$1")"
        PASS=$((PASS + 1))
    else
        echo "[FAIL] $1 not found in PATH"
        FAIL=$((FAIL + 1))
    fi
}

check_py_module() {
    if python3 -c "import $1" >/dev/null 2>&1; then
        echo "[OK]   python3 module: $1"
        PASS=$((PASS + 1))
    else
        echo "[FAIL] python3 module: $1 (not installed)"
        FAIL=$((FAIL + 1))
    fi
}

echo "================ Required commands ================"
check_cmd python3
check_cmd Rscript
check_cmd samtools
check_cmd xargs

echo "================ Python modules ================"
for m in numpy pandas scipy regex pyfaidx pysam xgboost sklearn; do
    check_py_module "$m"
done

echo "================ R ================"
if Rscript -e 'cat(as.character(getRversion()), "\n")' >/dev/null 2>&1; then
    echo "[OK]   Rscript runs"
    PASS=$((PASS + 1))
else
    echo "[FAIL] Rscript does not run"
    FAIL=$((FAIL + 1))
fi

echo "================ Reference genome ================"
REF="${MTDNAPIPE_REFERENCE:-$SCRIPT_DIR/reference/human_mtDNA.fasta}"
if [ -f "$REF" ]; then
    echo "[OK]   reference fasta: $REF"
    PASS=$((PASS + 1))
else
    echo "[FAIL] reference fasta not found: $REF"
    FAIL=$((FAIL + 1))
fi
if [ -f "$REF.fai" ]; then
    echo "[OK]   reference index (.fai) exists"
    PASS=$((PASS + 1))
else
    echo "[WARN] reference index (.fai) missing; run_pipeline.sh will create it with 'samtools faidx'"
fi

echo
echo "=================================================="
echo "Result: $PASS passed, $FAIL failed"
if [ "$FAIL" -gt 0 ]; then
    echo "Some requirements are missing. See README.md -> Installation."
    exit 1
fi
echo "Environment OK, you can run: bash run_pipeline.sh -i <sample_dir>"
