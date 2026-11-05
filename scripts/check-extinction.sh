#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."
repository_root="$PWD"

run_dir=$(mktemp -d "$PWD/build/extinction-check.XXXXXX")
trap 'rm -rf -- "$run_dir"' EXIT

(
  cd "$run_dir"
  "$repository_root/build/test-reversal-empty.exe" \
    --steps 5 --realisations 2 --start 3 --end 4 --seed 1
)

# Population includes time zero: six observations per realisation.
awk '
  NF != 6 { bad = 1 }
  {
    for (i = 1; i <= NF; i++)
      if ($i != 0) bad = 1
  }
  END {
    if (bad || NR != 2) {
      print "FAIL: empty-population output."
      exit 1
    }
  }
' "$run_dir/population.txt"

# Measurement window is steps 3, 4 and 5.
# No living-cell measurement exists, so all entries must be NaN.
awk '
  NF != 3 { bad = 1 }
  {
    for (i = 1; i <= NF; i++)
      if (tolower($i) !~ /^[-+]?nan(\([^)]*\))?$/)
        bad = 1
  }
  END {
    if (bad || NR != 2) {
      print "FAIL: extinction measurement window or missing values."
      exit 1
    }
  }
' "$run_dir/qav.txt"

awk '
  NF != 2 || $1 != NR || $2 != 0 { bad = 1 }
  END {
    if (bad || NR != 2) {
      print "FAIL: initially empty populations must have extinction time 0."
      exit 1
    }
  }
' "$run_dir/extinction.txt"

echo "PASS: empty initial population preserves output alignment."