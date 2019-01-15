#!/bin/sh

# Runs the benchmarks for various functions, prints results
for bench in *.g; do
    if [ "$bench" != "common.g" ]; then
        echo "Running $bench"
        gap -q <$bench
    fi
done
