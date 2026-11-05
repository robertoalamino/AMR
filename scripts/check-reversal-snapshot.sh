#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

test_dir=$(mktemp -d "$PWD/build/snapshot-check.XXXXXX")
trap 'rm -rf -- "$test_dir"' EXIT

printf '1\t-1\n2\t3\n3\t0\n' > "$test_dir/extinction.txt"

printf '0.1\t0.2\t0.3\n0.4\t0.5\t0.5\nnan\tnan\tnan\n' \
  > "$test_dir/qav.txt"

# Columns represent steps 2, 3 and 4.
bash scripts/extract-reversal-snapshot.sh \
  "$test_dir/qav.txt" "$test_dir/extinction.txt" 2 3 \
  > "$test_dir/snapshot.tsv"

awk '
  NR == 1 {
    if (NF != 6) bad = 1
    next
  }

  NF != 6 || $1 != NR - 1 || $2 != 3 { bad = 1 }

  NR == 2 {
    if ($3 != 0.2 || $4 != -1 || $5 != "alive" || $6 != 1)
      bad = 1
  }

  NR == 3 {
    if ($3 != 0.5 || $4 != 3 ||
        $5 != "extinction_step" || $6 != 1)
      bad = 1
  }

  NR == 4 {
    if (tolower($3) != "nan" || $4 != 0 ||
        $5 != "extinct_before" || $6 != 0)
      bad = 1
  }

  END {
    if (bad || NR != 4) {
      print "FAIL: snapshot values or extinction classifications."
      exit 1
    }
  }
' "$test_dir/snapshot.tsv"

# Step 5 is outside the supplied trajectory.
if bash scripts/extract-reversal-snapshot.sh \
    "$test_dir/qav.txt" "$test_dir/extinction.txt" 2 5 \
    >/dev/null 2>&1; then
  echo "FAIL: an out-of-range snapshot was accepted."
  exit 1
fi

echo "PASS: snapshot extraction, extinction status and missing values."