#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
bash scripts/build-historical.sh
# Compare exact content after normalising only Windows line endings.
compare_reference() {
  cmp <(sed 's/\r$//' "$1") <(sed 's/\r$//' "$2")
}
run_dir=$(mktemp -d "$PWD/build/historical-amdspace-check.XXXXXX")
trap 'rm -rf -- "$run_dir"' EXIT
(
  cd "$run_dir"
  ../historical/amdspace.exe --steps 30 --realisations 2 --seed 1
)
compare_reference tests/reference/amdspace-short/asymdrug.txt "$run_dir/asymdrug.txt"
echo "PASS: historical amdspace matches the saved reference exactly."
