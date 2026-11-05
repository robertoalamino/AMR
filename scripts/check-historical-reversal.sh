#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
bash scripts/build-historical.sh
# Compare exact content after normalising only Windows line endings.
compare_reference() {
  cmp <(sed 's/\r$//' "$1") <(sed 's/\r$//' "$2")
}
run_dir=$(mktemp -d "$PWD/build/historical-reversal-check.XXXXXX")
trap 'rm -rf -- "$run_dir"' EXIT
(
  cd "$run_dir"
  ../historical/reversaltimes.exe --steps 120 --realisations 2 --start 40 --end 60 --seed 1
)
compare_reference tests/reference/reversal-short/population.txt "$run_dir/population.txt"
compare_reference tests/reference/reversal-short/qav.txt "$run_dir/qav.txt"
echo "PASS: historical reversal matches the saved reference exactly."
