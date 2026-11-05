#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

check_valid() {
  local expected="$1"
  shift
  local actual

  actual=$(./build/test-config.exe "$@")

  if [[ "$actual" != "$expected" ]]; then
    echo "FAIL: seed arguments '$*': expected $expected, got $actual."
    return 1
  fi
}

check_invalid() {
  local status=0

  ./build/test-config.exe "$@" >/dev/null 2>&1 || status=$?

  if [[ "$status" -ne 1 ]]; then
    echo "FAIL: invalid arguments '$*' returned $status; expected 1."
    return 1
  fi
}

check_valid "-1 120 2 40 60"
check_valid "-1 120 2 40 60" --seed 1
check_valid "-42 120 2 40 60" --seed 42

check_valid "-7 80 3 10 30" \
  --steps 80 --realisations 3 --start 10 --end 30 --seed 7

# End is exclusive: treatment can continue through the final step.
check_valid "-1 120 2 40 121" --end 121

# Validation must not depend on argument order.
check_valid "-1 20 2 5 15" --end 15 --start 5 --steps 20

check_invalid --seed
check_invalid --seed 0
check_invalid --seed -1
check_invalid --seed abc
check_invalid --seed 1.5
check_invalid --seed 42junk
check_invalid --seed 999999999999999999999999999999999999
check_invalid --unknown 42
check_invalid --seed 42 extra
check_invalid --steps 0
check_invalid --steps -1
check_invalid --steps abc
check_invalid --steps 999999999999999999999999
check_invalid --realisations 0
check_invalid --start 0
check_invalid --start 60 --end 60
check_invalid --start 70 --end 60
check_invalid --end 122
check_invalid --steps 20

echo "PASS: runtime settings and invalid arguments."