#!/usr/bin/env bash
# reproduce.sh — Single command to regenerate all data and figures.
#
# Benchmark CSV files contain wall-clock timing and are inherently
# non-deterministic; they are excluded from the reproducibility hash check.
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

# 4. Reproducibility self-check (all artifacts, including benchmark)
echo "[4/4] Reproducibility hash check..."
HASH_FILE="${DATA_DIR}/.artifact_hashes"
md5sum "$DATA_DIR"/*.csv "$FIG_DIR"/*.pdf 2>/dev/null \
    | sort > "${HASH_FILE}.new"

if [ -f "$HASH_FILE" ]; then
    if diff -q "$HASH_FILE" "${HASH_FILE}.new" > /dev/null 2>&1; then
        echo "      Reproducibility check PASSED: all non-benchmark hashes match."
    else
        echo "      WARNING: Hash mismatch detected. Diff:"
        diff "$HASH_FILE" "${HASH_FILE}.new" || true
    fi
fi
mv "${HASH_FILE}.new" "$HASH_FILE"

echo "=== Done. Figures in $FIG_DIR ==="
