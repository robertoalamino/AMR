#!/usr/bin/env bash
set -euo pipefail

# Run from the repository root, regardless of where this script is launched.
cd "$(dirname "$0")/.."

mkdir -p build

echo "Building antimicrobial parameter sweep..."
gcc -std=c99 -Wall -Wextra -Wpedantic -O2 \
  src/amdspace.c src/genotype.c src/random.c src/response.c src/lattice.c src/config.c src/population.c \
  -o build/amdspace.exe -lm

echo "Building reversal simulation with original default settings..."
gcc -std=c99 -Wall -Wextra -Wpedantic -O2 \
  src/reversaltimes.c src/genotype.c src/random.c src/response.c src/lattice.c src/config.c src/population.c \
  -o build/reversaltimes.exe -lm

echo "Building lattice tests..."
gcc -std=c99 -Wall -Wextra -Wpedantic -O2 -Isrc \
  tests/test_lattice.c src/lattice.c -o build/test-lattice.exe

echo "Building seed parser test..."
gcc -std=c99 -Wall -Wextra -Wpedantic -O2 -Isrc \
  tests/test_config.c src/config.c -o build/test-config.exe

echo "Building empty neighbour tests..."
gcc -std=c99 -Wall -Wextra -Wpedantic -O2 -Isrc \
  tests/test_empty_neighbours.c src/lattice.c \
  -o build/test-empty-neighbours.exe

echo "Building population tests..."
gcc -std=c99 -Wall -Wextra -Wpedantic -O2 -Isrc \
  tests/test_population.c src/population.c \
  -o build/test-population.exe

echo "Building empty-population reversal test..."
gcc -std=c99 -Wall -Wextra -Wpedantic -O2 -Isrc \
  src/reversaltimes.c src/genotype.c src/response.c \
  src/lattice.c src/config.c src/population.c \
  tests/empty_population_random.c \
  -o build/test-reversal-empty.exe -lm

echo "Building empty-population sweep test..."
gcc -std=c99 -Wall -Wextra -Wpedantic -O2 -Isrc \
  src/amdspace.c src/genotype.c src/response.c \
  src/lattice.c src/config.c src/population.c \
  tests/empty_population_random.c \
  -o build/test-amdspace-empty.exe -lm

echo "Build completed."