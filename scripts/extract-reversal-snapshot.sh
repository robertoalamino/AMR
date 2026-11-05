#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 4 ||
      ! "$3" =~ ^[1-9][0-9]*$ ||
      ! "$4" =~ ^[1-9][0-9]*$ ]]; then
  echo "Usage: bash scripts/extract-reversal-snapshot.sh QAV_FILE EXTINCTION_FILE START_STEP TARGET_STEP" >&2
  exit 1
fi

awk -v start="$3" -v target="$4" '
  BEGIN {
    OFS = "\t"
    column = target - start + 1

    if (column < 1) {
      print "Target step precedes the measurement window." > "/dev/stderr"
      bad = 1
      exit 1
    }
  }

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

  FILENAME == ARGV[1] {
    if (NF != 2 || $1 != FNR ||
        $2 !~ /^-?[0-9]+$/ || $2 < -1) {
      fail("invalid extinction record on row " FNR)
      next
    }

    extinction[FNR] = $2 + 0
    extinction_count++
    next
  }

  FILENAME == ARGV[2] {
    trajectory_count++

    if (FNR == 1)
      columns = NF

    if (NF != columns || column > NF || FNR > extinction_count) {
      fail("trajectory dimensions do not match on row " FNR)
      next
    }

    value = $column

    if (!is_nan(value) &&
        (!is_number(value) || value < 0 || value > 1)) {
      fail("invalid susceptibility on row " FNR)
      next
    }

    values[FNR] = value
  }

  END {
    if (bad)
      exit 1

    if (trajectory_count == 0 ||
        trajectory_count != extinction_count) {
      print "FAIL: missing or mismatched realisations." > "/dev/stderr"
      exit 1
    }

    for (row = 1; row <= trajectory_count; row++) {
      if (extinction[row] > start + columns - 1) {
        print "FAIL: extinction time exceeds the run." > "/dev/stderr"
        exit 1
      }
    }

    print "realisation", "step", "susceptibility", \
          "extinction_time", "status", "available"

    for (row = 1; row <= trajectory_count; row++) {
      time = extinction[row]

      if (time == -1 || time > target)
        status = "alive"
      else if (time == target)
        status = "extinction_step"
      else
        status = "extinct_before"

      print row, target, values[row], time, status, \
            is_nan(values[row]) ? 0 : 1
    }
  }
' "$2" "$1"