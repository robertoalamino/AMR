#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p build/historical
for program in amdspace reversaltimes; do
  echo "Building historical $program..."
  gcc -std=c99 -Wall -Wextra -Wpedantic -O2 \
    "historical/src/$program.c" historical/src/genotype.c \
    historical/src/random.c historical/src/response.c \
    historical/src/lattice.c historical/src/config.c \
    -o "build/historical/$program.exe" -lm
done
