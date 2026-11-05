#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."
repository_root="$PWD"

run_dir=$(mktemp -d "$PWD/build/amdspace-extinction.XXXXXX")
trap 'rm -rf -- "$run_dir"' EXIT

(
  cd "$run_dir"
  "$repository_root/build/test-amdspace-empty.exe" \
    --steps 3 --realisations 2 --seed 1
)

# Each parameter combination has two initially empty realisations.
awk '
  {
    expected_realisation = (NR - 1) % 2 + 1

    if (NF != 4 || $3 != expected_realisation || $4 != 0)
      bad = 1
  }

  END {
    if (bad || NR != 242) {
      print "FAIL: expected 242 initially empty realisations."
      exit 1
    }
  }
' "$run_dir/extinction.txt"

# No valid death-probability measurement exists in any realisation.
awk '
  NF != 3 || tolower($3) !~ /^[-+]?nan(\([^)]*\))?$/ { bad = 1 }

  END {
    if (bad || NR != 121) {
      print "FAIL: expected 121 unavailable sweep measurements."
      exit 1
    }
  }
' "$run_dir/asymdrug.txt"

bash scripts/check-sweep-records.sh "$run_dir" 2 3

echo "PASS: sweep records initial extinction and unavailable measurements."