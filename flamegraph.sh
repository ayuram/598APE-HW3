#!/bin/bash

set -euo pipefail

# === CONFIGURATION ===
EXECUTABLE="./main.exe"
NPLANETS=1000
TIMESTEPS=5000
PERF_DATA="perf.data"
PERF_SCRIPT="perf_unfolded.perf"
FLAMEGRAPH_SVG="flamegraph.svg"
SUMMARY_TXT="perf_summary.txt"

# === CHECK DEPENDENCIES ===
for cmd in make perf git gcc; do
    command -v "$cmd" >/dev/null 2>&1 || { echo "$cmd is not installed. Aborting." >&2; exit 1; }
done

# === BUILD PROJECT ===
echo "Building the project..."
make -j

# === PERF RECORD ===
echo "Running perf record on $EXECUTABLE..."
perf record -F 99 -g -- "$EXECUTABLE" "$NPLANETS" "$TIMESTEPS"

# === PERF REPORT (Interactive) ===
echo "Generating interactive perf report..."
perf report

# === GENERATE SUMMARY TEXT REPORT ===
echo "Generating summary report (CPU cycles, instructions, cache misses, etc.)..."

perf stat -d -- "$EXECUTABLE" "$NPLANETS" "$TIMESTEPS" | tee "$SUMMARY_TXT"

echo "Summary report saved to $SUMMARY_TXT 📄"

# === FLAMEGRAPH ===
echo "Generating flamegraph..."
if [ ! -d "FlameGraph" ]; then
    echo "FlameGraph repo not found, cloning..."
    git clone https://github.com/brendangregg/FlameGraph.git
fi

# Expand the raw perf data
perf script > "$PERF_SCRIPT"

# Collapse the stacks
./FlameGraph/stackcollapse-perf.pl "$PERF_SCRIPT" > out.folded

# Generate the flamegraph
./FlameGraph/flamegraph.pl out.folded > "$FLAMEGRAPH_SVG"

echo "Flamegraph generated at $FLAMEGRAPH_SVG ✅"

# === FINAL NOTE ===
if command -v xdg-open &> /dev/null; then
    echo "Opening flamegraph in default browser..."
    xdg-open "$FLAMEGRAPH_SVG"
elif command -v open &> /dev/null; then  # macOS
    echo "Opening flamegraph in default browser..."
    open "$FLAMEGRAPH_SVG"
else
    echo "Please manually open $FLAMEGRAPH_SVG in your browser 🌍"
fi

echo "All done! 🎯"
