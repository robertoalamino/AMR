#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

summary_file=$(mktemp "$PWD/build/sweep-summary.XXXXXX")
trap 'rm -f -- "$summary_file"' EXIT

# Three parameter pairs:
# 1. Mixed alive, extinct and unavailable.
# 2. All extinct, with valid last-observed values.
# 3. All unavailable.
printf '%s\n' \
  '0.0 0.0 1 0.2 -1' \
  '0.0 0.0 2 0.4 3' \
  '0.0 0.0 3 nan 0' \
  '0.0 0.2 1 0.2 3' \
  '0.0 0.2 2 0.6 4' \
  '0.0 0.4 1 nan 0' \
  '0.0 0.4 2 nan 0' |
  bash scripts/summarise-sweep.sh > "$summary_file"

awk '
  function approximately_equal(a, b) {
    return (a - b < 0.000001) && (b - a < 0.000001)
  }

  NR == 2 {
    if ($3 != 3 || $4 != 1 || $5 != 2 || $6 != 1 ||
        $7 != "NaN" || $8 != 2 || !approximately_equal($9, 0.3) ||
        $10 != 1 || !approximately_equal($11, 0.2))
      bad = 1
  }

  NR == 3 {
    if ($3 != 2 || $4 != 0 || $5 != 2 || $6 != 0 ||
        !approximately_equal($7, 0.4) || $8 != 2 || !approximately_equal($9, 0.4) ||
        $10 != 0 || $11 != "NaN")
      bad = 1
  }

  NR == 4 {
    if ($3 != 2 || $4 != 0 || $5 != 2 || $6 != 2 ||
        $7 != "NaN" || $8 != 0 || $9 != "NaN" ||
        $10 != 0 || $11 != "NaN")
      bad = 1
  }

  END {
    if (bad || NR != 4) {
      print "FAIL: sweep summary values."
      exit 1
    }
  }
' "$summary_file"

# Invalid probabilities must be rejected.
if printf '0.0 0.0 1 1.5 -1\n' |
   bash scripts/summarise-sweep.sh >/dev/null 2>&1; then
  echo "FAIL: invalid susceptibility was accepted."
  exit 1
fi

echo "PASS: sweep summaries and invalid input."