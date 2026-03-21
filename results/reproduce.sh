#!/usr/bin/env bash
# reproduce.sh — Single command to regenerate all data and figures.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR="${SCRIPT_DIR}/build"
DATA_DIR="${SCRIPT_DIR}/data"
FIG_DIR="${SCRIPT_DIR}/figures"

echo "=== Results Reproduction Pipeline ==="

# 1. Configure and build the data generator
echo "[1/4] Building data generator..."
cmake -S "$SCRIPT_DIR" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release 2>&1 | tail -3
cmake --build "$BUILD_DIR" --parallel 2>&1 | tail -3
echo "      Build complete."

# 2. Run the data generator
echo "[2/4] Generating experiment data..."
mkdir -p "$DATA_DIR" "$FIG_DIR"
cd "$BUILD_DIR"
LD_LIBRARY_PATH="${SCRIPT_DIR}/../build/ql:${LD_LIBRARY_PATH:-}" ./generate_data "$DATA_DIR"
cd "$SCRIPT_DIR"

# 3. Generate figures
echo "[3/4] Generating figures..."
python3 "$SCRIPT_DIR/plot_figures.py"

# 4. Reproducibility self-check (split: deterministic vs non-deterministic)
echo "[4/4] Reproducibility hash check..."
HASH_FILE="${DATA_DIR}/.artifact_hashes"

# Deterministic artifacts only (exclude benchmark CSVs/PDFs with wall-clock timing)
{
  find "$DATA_DIR" -name '*.csv' ! -name 'benchmark_*.csv' -print0 | xargs -0 md5sum 2>/dev/null
  find "$FIG_DIR" -name '*.pdf' ! -name 'fig8_benchmark.pdf' -print0 | xargs -0 md5sum 2>/dev/null
} | sort > "${HASH_FILE}.new"

if [ -f "$HASH_FILE" ]; then
    if diff -q "$HASH_FILE" "${HASH_FILE}.new" > /dev/null 2>&1; then
        echo "      Reproducibility check PASSED: all deterministic hashes match."
    else
        echo "      WARNING: Hash mismatch detected. Diff:"
        diff "$HASH_FILE" "${HASH_FILE}.new" || true
    fi
fi
mv "${HASH_FILE}.new" "$HASH_FILE"

echo "=== Done. Figures in $FIG_DIR ==="
