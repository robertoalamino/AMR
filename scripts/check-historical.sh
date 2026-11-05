#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
bash scripts/build-historical.sh
# Compare exact content after normalising only Windows line endings.
compare_reference() {
  cmp <(sed 's/\r$//' "$1") <(sed 's/\r$//' "$2")
}
run_dir=$(mktemp -d "$PWD/build/historical-check.XXXXXX")
trap 'rm -rf -- "$run_dir"' EXIT
mkdir "$run_dir/reversal" "$run_dir/sweep"
(
 cd "$run_dir/reversal"
 ../../historical/reversaltimes.exe --steps 120 --realisations 2 --start 40 --end 60 --seed 1
)
for file in population.txt qav.txt; do
 compare_reference "tests/reference/reversal-short/$file" "$run_dir/reversal/$file"
done
(
 cd "$run_dir/sweep"
 ../../historical/amdspace.exe --steps 30 --realisations 2 --seed 1
)
compare_reference tests/reference/amdspace-short/asymdrug.txt "$run_dir/sweep/asymdrug.txt"
echo "PASS: both historical regressions match saved references exactly."
