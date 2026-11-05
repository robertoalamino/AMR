#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

bash scripts/build.sh

./build/test-lattice.exe
./build/test-empty-neighbours.exe
./build/test-population.exe

bash scripts/check-config.sh
bash scripts/check-output-errors.sh
bash scripts/check-simulations.sh
bash scripts/check-extinction.sh
bash scripts/check-extinction-during-run.sh
bash scripts/check-amdspace-extinction.sh
bash scripts/check-amdspace-extinction-during-run.sh
bash scripts/check-sweep-summary.sh
bash scripts/check-reversal-summary.sh
bash scripts/check-reversal-snapshot.sh
bash scripts/check-reversal-bins.sh
bash scripts/check-parent-selection.sh

echo "PASS: all corrected-implementation checks completed."