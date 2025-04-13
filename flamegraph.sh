#!/bin/bash

set -e  # Exit immediately if a command fails

TARGET=main.exe
PERF_REPORT=perf_report.txt

# Compile the program
echo "Compiling..."
g++ -O3 -g -fno-omit-frame-pointer -Werror -fopenmp main.cpp -o $TARGET -lm

# Run the program in the background
echo "Running simulation in background..."
./$TARGET 1000 5000 &
PID=$!

# Wait a bit for it to start up
sleep 1

# Collect perf stat (per-thread)
echo "Collecting perf stat..."
perf stat --per-thread -p $PID -o $PERF_REPORT

# Wait for program to finish
wait $PID

# Now do a perf record + flamegraph
echo "Generating flamegraph..."
perf record -F 99 -g ./$TARGET 1000 5000
perf script | ~/Flamegraph/stackcollapse-perf.pl | ~/Flamegraph/flamegraph.pl > flamegraph.svg

echo "Done!"
