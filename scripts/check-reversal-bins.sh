#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
mkdir -p "$repo_dir/build"
test_dir=$(mktemp -d "$repo_dir/build/reversal-bins.XXXXXX")
trap 'rm -rf "$test_dir"' EXIT

cat > "$test_dir/input.tsv" <<'EOF'
realisation step susceptibility extinction_time status available
1 5 0 -1 alive 1
2 5 0.01 -1 alive 1
3 5 0.29 3 extinct_before 1
4 5 0.99 5 extinction_step 1
5 5 1 -1 alive 1
6 5 nan 0 extinct_before 0
EOF

bash "$repo_dir/scripts/bin-reversal-snapshot.sh" \
  < "$test_dir/input.tsv" \
  > "$test_dir/bins.tsv" \
  2> "$test_dir/counts.txt"

awk '
  function approximately_equal(a, b) {
    return (a-b < 1e-12 && b-a < 1e-12)
  }

  NR == 1 {
    if (NF != 10 || $1 != "lower" || $10 != "fraction_alive_available")
      bad = 1
    next
  }

  {
    bin = NR - 2

    expected_available = (bin == 0 || bin == 1 || bin == 29) ? 1 : 0
    if (bin == 99) expected_available = 2

    expected_alive = (bin == 0 || bin == 1 || bin == 99) ? 1 : 0
    expected_extinct = (bin == 29) ? 1 : 0
    expected_at_step = (bin == 99) ? 1 : 0

    if (NF != 10 ||
        !approximately_equal($1, bin/100) ||
        !approximately_equal($2, (bin+1)/100) ||
        !approximately_equal($3, (bin+0.5)/100) ||
        $4 != expected_available ||
        $5 != expected_alive ||
        $6 != expected_extinct ||
        $7 != expected_at_step ||
        !approximately_equal($8, expected_available/6) ||
        !approximately_equal($9, expected_available/5) ||
        !approximately_equal($10, expected_alive/3)) {
      print "FAIL: histogram bin", bin
      bad = 1
    }

    total_fraction += $8
    available_fraction += $9
    alive_fraction += $10
  }

  END {
    if (NR != 101 ||
        !approximately_equal(total_fraction, 5/6) ||
        !approximately_equal(available_fraction, 1) ||
        !approximately_equal(alive_fraction, 1))
      bad = 1

    if (bad) {
      print "FAIL: histogram boundaries, categories or normalisation."
      exit 1
    }
  }
' "$test_dir/bins.tsv"

# With no available measurements, conditional fractions are undefined.
printf 'realisation step susceptibility extinction_time status available\n1 5 nan 0 extinct_before 0\n' |
  bash "$repo_dir/scripts/bin-reversal-snapshot.sh" \
    > "$test_dir/unavailable.tsv" \
    2> "$test_dir/unavailable-counts.txt"

awk '
  NR == 1 { next }
  {
    if ($4 != 0 || $8 != 0 ||
        tolower($9) != "nan" || tolower($10) != "nan")
      bad = 1
  }
  END {
    if (NR != 101 || bad) {
      print "FAIL: unavailable-only histogram."
      exit 1
    }
  }
' "$test_dir/unavailable.tsv"

# Measurements outside [0, 1] must be rejected.
if printf 'realisation step susceptibility extinction_time status available\n1 5 1.1 -1 alive 1\n' |
  bash "$repo_dir/scripts/bin-reversal-snapshot.sh" \
    > "$test_dir/invalid.tsv" 2> "$test_dir/invalid.log"
then
  echo "FAIL: out-of-range susceptibility was accepted."
  exit 1
fi

echo "PASS: reversal histogram boundaries, categories and normalisation."