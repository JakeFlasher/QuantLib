#!/usr/bin/env bash
# reproduce.sh — Single command to regenerate all data and figures.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR="${SCRIPT_DIR}/build"
DATA_DIR="${SCRIPT_DIR}/data"
FIG_DIR="${SCRIPT_DIR}/figures"

echo "=== Results Reproduction Pipeline ==="

# 1. Configure and build the data generator
echo "[1/3] Building data generator..."
cmake -S "$SCRIPT_DIR" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release 2>&1 | tail -3
cmake --build "$BUILD_DIR" --parallel 2>&1 | tail -3
echo "      Build complete."

# 2. Run the data generator
echo "[2/3] Generating experiment data..."
mkdir -p "$DATA_DIR" "$FIG_DIR"
cd "$DATA_DIR"
"$BUILD_DIR/generate_data"
cd "$SCRIPT_DIR"

# 3. Generate figures
echo "[3/3] Generating figures..."
python3 "$SCRIPT_DIR/plot_figures.py"

echo "=== Done. Figures in $FIG_DIR ==="
