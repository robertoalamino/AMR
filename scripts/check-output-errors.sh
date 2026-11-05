#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

test_dir=$(mktemp -d "$PWD/build/output-errors.XXXXXX")
trap 'rm -rf -- "$test_dir"' EXIT

check_failure() {
  local name="$1"
  local executable="$2"
  local blocked_file="$3"
  local expected_message="$4"
  local case_dir="$test_dir/$name"
  local status=0

  mkdir -p "$case_dir/$blocked_file"

  (
    cd "$case_dir"
    "$executable"
  ) >"$case_dir/stdout.txt" 2>"$case_dir/stderr.txt" || status=$?

  if [[ "$status" -ne 1 ]]; then
    echo "FAIL: $name returned $status; expected 1."
    return 1
  fi

  if ! grep -Fq "$expected_message" "$case_dir/stderr.txt"; then
    echo "FAIL: $name did not report the expected error."
    cat "$case_dir/stderr.txt"
    return 1
  fi

  echo "PASS: $name."
}

check_failure \
  "sweep-output" \
  "$PWD/build/amdspace.exe" \
  "asymdrug.txt" \
  "Cannot open asymdrug.txt"

check_failure \
  "reversal-population" \
  "$PWD/build/reversaltimes.exe" \
  "population.txt" \
  "Cannot open population.txt"

check_failure \
  "reversal-probability" \
  "$PWD/build/reversaltimes.exe" \
  "qav.txt" \
  "Cannot open qav.txt"

check_failure \
  "reversal-extinction" \
  "$PWD/build/reversaltimes.exe" \
  "extinction.txt" \
  "Cannot open extinction.txt"

check_failure \
  "sweep-extinction" \
  "$PWD/build/amdspace.exe" \
  "extinction.txt" \
  "Cannot open extinction.txt"

check_failure \
  "sweep-realisations" \
  "$PWD/build/amdspace.exe" \
  "sweep-realisations.txt" \
  "Cannot open sweep-realisations.txt"