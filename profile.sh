#!/bin/bash

set -e  # Exit immediately on error

# Config
TARGET="main.exe"
FLAMEGRAPH_DIR=~/Flamegraph
ARGS="1000 5000"  # Arguments for your n-body program

# Step 1: Build
echo "🛠️  Building the project..."
make

# Step 2: Clean old perf artifacts
echo "🧹 Cleaning old perf data..."
rm -f perf.data perf.data.old perf_report.txt perf_stat.txt flamegraph.svg

# Step 4: Collect CPU cycles and branch misses
echo "🧠 Collecting CPU and branch stats (perf_stat.txt)..."
perf stat -e cycles,instructions,branches,branch-misses,cache-references,cache-misses -o perf_stat.txt ./$TARGET $ARGS

echo "✅ Done! Files generated:"
echo " - perf_stat.txt (cycles, branch misses, cache stats)"
