#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: bash scripts/check-sweep-records.sh DIRECTORY REALISATIONS STEPS" >&2
  exit 1
fi

run_dir="$1"
realisations="$2"
steps="$3"

awk -v realisations="$realisations" -v steps="$steps" '
  function is_nan(value) {
    return tolower(value) ~ /^[-+]?nan(\([^)]*\))?$/
  }

  function is_number(value) {
    return value ~ /^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/
  }

  function fail(message) {
    print "FAIL: " message > "/dev/stderr"
    bad = 1
  }

  # First file: one aggregate for each parameter combination.
  FILENAME == ARGV[1] {
    if (NF != 3) {
      fail("invalid aggregate row " FNR)
      next
    }

    alpha[FNR] = $1
    beta[FNR] = $2
    aggregate[FNR] = $3
    aggregate_count++
    next
  }

  # Second file: one extinction record for each realisation.
  FILENAME == ARGV[2] {
    if (NF != 4) {
      fail("invalid extinction row " FNR)
      next
    }

    ext_alpha[FNR] = $1
    ext_beta[FNR] = $2
    ext_number[FNR] = $3
    ext_time[FNR] = $4
    extinction_count++
    next
  }

  # Third file: individual susceptibility records.
  FILENAME == ARGV[3] {
    record_count++

    if (NF != 5) {
      fail("invalid susceptibility row " FNR)
      next
    }

    pair = int((FNR - 1) / realisations) + 1
    expected_number = (FNR - 1) % realisations + 1

    if (pair > aggregate_count ||
        $1 != alpha[pair] || $2 != beta[pair] ||
        $1 != ext_alpha[FNR] || $2 != ext_beta[FNR] ||
        $3 != expected_number || $3 != ext_number[FNR] ||
        $5 != ext_time[FNR]) {
      fail("parameter or realisation mismatch on row " FNR)
      next
    }

    if ($5 !~ /^-?[0-9]+$/ || $5 < -1 || $5 > steps) {
      fail("invalid extinction time on row " FNR)
      next
    }

    if (is_nan($4)) {
      missing[pair] = 1
    } else if (!is_number($4) || $4 < 0 || $4 > 1) {
      fail("invalid susceptibility on row " FNR)
    } else {
      total[pair] += $4
    }
  }

  END {
    expected_records = 121 * realisations

    if (aggregate_count != 121 ||
        extinction_count != expected_records ||
        record_count != expected_records)
      fail("unexpected number of rows")

    for (pair = 1; pair <= aggregate_count; pair++) {
      if (missing[pair]) {
        # The current C aggregate propagates NaN.
        if (!is_nan(aggregate[pair]))
          fail("expected NaN aggregate for parameter pair " pair)
      } else if (!is_number(aggregate[pair])) {
        fail("invalid aggregate for parameter pair " pair)
      } else {
        expected = total[pair] / realisations
        difference = aggregate[pair] - expected
        if (difference < 0) difference = -difference

        # asymdrug.txt is printed to six decimal places.
        if (difference > 0.00000051)
          fail("aggregate mismatch for parameter pair " pair)
      }
    }

    if (bad)
      exit 1

    print "PASS: sweep realisation records and aggregate reconstruction."
  }
' "$run_dir/asymdrug.txt" \
  "$run_dir/extinction.txt" \
  "$run_dir/sweep-realisations.txt"