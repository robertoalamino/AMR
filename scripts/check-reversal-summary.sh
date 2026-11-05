#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

test_dir=$(mktemp -d "$PWD/build/reversal-summary-check.XXXXXX")
trap 'rm -rf -- "$test_dir"' EXIT

printf '1\t-1\n2\t3\n3\t0\n' > "$test_dir/extinction.txt"

printf '0.1\t0.2\t0.3\n0.4\t0.5\t0.5\nnan\tnan\tnan\n' \
  > "$test_dir/qav.txt"

bash scripts/summarise-reversal.sh \
  "$test_dir/qav.txt" "$test_dir/extinction.txt" 2 \
  > "$test_dir/summary.tsv"

awk '
  function approximately_equal(a, b) {
    return (a - b < 0.000001) && (b - a < 0.000001)
  }

  NR == 2 {
    if ($1 != 2 || $2 != 3 || $3 != 2 || $4 != 1 ||
        $5 != 0 || $6 != 1 || $7 != "NaN" ||
        $8 != 2 || !approximately_equal($9, 0.25) ||
        $10 != 2 || !approximately_equal($11, 0.25) ||
        $12 != 0)
      bad = 1
  }

  NR == 3 {
    if ($1 != 3 || $2 != 3 || $3 != 1 || $4 != 1 ||
        $5 != 1 || $6 != 1 || $7 != "NaN" ||
        $8 != 2 || !approximately_equal($9, 0.35) ||
        $10 != 1 || !approximately_equal($11, 0.2) ||
        $12 != 0)
      bad = 1
  }

  NR == 4 {
    if ($1 != 4 || $2 != 3 || $3 != 1 || $4 != 2 ||
        $5 != 0 || $6 != 1 || $7 != "NaN" ||
        $8 != 2 || !approximately_equal($9, 0.4) ||
        $10 != 1 || !approximately_equal($11, 0.3) ||
        $12 != 1)
      bad = 1
  }

  END {
    if (bad || NR != 4) {
      print "FAIL: reversal summary values."
      exit 1
    }
  }
' "$test_dir/summary.tsv"

# A missing extinction record must cause an error.
printf '1\t-1\n2\t3\n' > "$test_dir/extinction-incomplete.txt"

if bash scripts/summarise-reversal.sh \
    "$test_dir/qav.txt" "$test_dir/extinction-incomplete.txt" 2 \
    >/dev/null 2>&1; then
  echo "FAIL: mismatched realisation counts were accepted."
  exit 1
fi

echo "PASS: reversal summaries and mismatched input."