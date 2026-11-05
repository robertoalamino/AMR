#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

test_dir=$(mktemp -d "$PWD/build/simulation-check.XXXXXX")
trap 'rm -rf -- "$test_dir"' EXIT

check_table() {
  local file="$1"
  local rows="$2"
  local columns="$3"
  local kind="$4"

  awk -v rows="$rows" -v columns="$columns" -v kind="$kind" '
    function fail(message) {
      print "FAIL: " FILENAME ": " message > "/dev/stderr"
      bad = 1
      exit 1
    }

    {
      if (NF != columns)
        fail("unexpected column count on row " NR)

      for (i = 1; i <= NF; i++) {
        if ($i !~ /^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/)
          fail("non-finite or invalid value on row " NR)

        value = $i + 0

        if (kind == "population" &&
            (value < 0 || value > 2500 || value != int(value)))
          fail("invalid population on row " NR)

        if (kind == "probability" && (value < 0 || value > 1))
          fail("invalid probability on row " NR)

        if (kind == "sweep") {
          if (i <= 2 && (value < -1 || value > 1))
            fail("invalid antimicrobial parameter on row " NR)
          if (i == 3 && (value < 0 || value > 1))
            fail("invalid probability on row " NR)
        }
      }
    }

    END {
      if (!bad && NR != rows) {
        print "FAIL: " FILENAME ": unexpected row count" > "/dev/stderr"
        exit 1
      }
    }
  ' "$file"
}

mkdir -p "$test_dir/reversal-1" "$test_dir/reversal-2"
mkdir -p "$test_dir/sweep-1" "$test_dir/sweep-2"

for repeat in 1 2; do
  (
    cd "$test_dir/reversal-$repeat"
    "$OLDPWD/build/reversaltimes.exe" \
      --steps 120 --realisations 2 --start 40 --end 60 --seed 1
  )

  (
    cd "$test_dir/sweep-$repeat"
    "$OLDPWD/build/amdspace.exe" \
      --steps 30 --realisations 2 --seed 1
  )
done

check_table "$test_dir/reversal-1/population.txt" 2 121 population
check_table "$test_dir/reversal-1/qav.txt" 2 81 probability
check_table "$test_dir/sweep-1/asymdrug.txt" 121 3 sweep

cmp "$test_dir/reversal-1/population.txt" \
    "$test_dir/reversal-2/population.txt"

cmp "$test_dir/reversal-1/qav.txt" \
    "$test_dir/reversal-2/qav.txt"

cmp "$test_dir/sweep-1/asymdrug.txt" \
    "$test_dir/sweep-2/asymdrug.txt"

cmp "$test_dir/reversal-1/extinction.txt" \
    "$test_dir/reversal-2/extinction.txt"

cmp "$test_dir/sweep-1/extinction.txt" \
    "$test_dir/sweep-2/extinction.txt"

# Each sweep row has two corresponding realisation records.
awk '
  FNR == NR {
    alpha[NR] = $1
    beta[NR] = $2
    combinations = NR
    next
  }

  {
    record = FNR
    pair = int((record - 1) / 2) + 1
    expected_realisation = (record - 1) % 2 + 1

    if (NF != 4 ||
        pair > combinations ||
        $1 != alpha[pair] || $2 != beta[pair] ||
        $3 != expected_realisation ||
        $4 !~ /^-?[0-9]+$/ ||
        $4 < -1 || $4 > 30)
      bad = 1
  }

  END {
    if (bad || combinations != 121 || FNR != 242) {
      print "FAIL: invalid sweep extinction records."
      exit 1
    }
  }
' "$test_dir/sweep-1/asymdrug.txt" \
  "$test_dir/sweep-1/extinction.txt"

cmp "$test_dir/sweep-1/sweep-realisations.txt" \
    "$test_dir/sweep-2/sweep-realisations.txt"

bash scripts/check-sweep-records.sh "$test_dir/sweep-1" 2 30

echo "PASS: simulation output structure, ranges and repeatability."