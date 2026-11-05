#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 || ! "$3" =~ ^[1-9][0-9]*$ ]]; then
  echo "Usage: bash scripts/summarise-reversal.sh QAV_FILE EXTINCTION_FILE START_STEP" >&2
  exit 1
fi

qav_file="$1"
extinction_file="$2"
start_step="$3"

awk -v start="$start_step" '
  BEGIN {
    OFS = "\t"
  }

  function is_nan(value) {
    return tolower(value) ~ /^[-+]?nan(\([^)]*\))?$/
  }

  function is_number(value) {
    return value ~ /^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/
  }

  function mean(sum, count) {
    return count ? sprintf("%.17g", sum / count) : "NaN"
  }

  function fail(message) {
    print "FAIL: " message > "/dev/stderr"
    bad = 1
  }

  # First file: realisation number and extinction time.
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

  # Second file: one susceptibility trajectory per realisation.
  FILENAME == ARGV[2] {
    row = FNR
    trajectory_count++

    if (row > extinction_count) {
      fail("more trajectories than extinction records")
      next
    }

    if (row == 1)
      columns = NF
    else if (NF != columns) {
      fail("inconsistent trajectory length on row " row)
      next
    }

    for (i = 1; i <= NF; i++) {
      step = start + i - 1
      time = extinction[row]

      total[i]++

      if (time == -1 || time > step)
        alive[i]++
      else if (time == step)
        extinct_this_step[i]++
      else
        extinct_before[i]++

      if (is_nan($i)) {
        unavailable[i]++
      } else if (!is_number($i) || $i < 0 || $i > 1) {
        fail("invalid susceptibility in row " row ", column " i)
      } else {
        available_count[i]++
        available_sum[i] += $i

        if (time == -1 || time > step) {
          alive_valid[i]++
          alive_sum[i] += $i
        }

        if (time >= 0 && time < step)
          carried_count[i]++
      }
    }
  }

  END {
    if (extinction_count == 0 ||
        trajectory_count != extinction_count || columns == 0)
      fail("missing or mismatched realisations")

    final_step = start + columns - 1

    for (row = 1; row <= extinction_count; row++)
      if (extinction[row] > final_step)
        fail("extinction time exceeds the observation period")

    if (bad)
      exit 1

    print "step", "total", "alive_after_step", \
          "extinct_before_step", "extinction_at_step", \
          "unavailable", "historical_style_mean", \
          "available_count", "available_mean", \
          "alive_valid_count", "alive_mean", \
          "carried_after_extinction_count"

    for (i = 1; i <= columns; i++) {
      historical_mean = unavailable[i] \
        ? "NaN" : mean(available_sum[i], available_count[i])

      print start + i - 1, total[i], alive[i] + 0, \
            extinct_before[i] + 0, extinct_this_step[i] + 0, \
            unavailable[i] + 0, historical_mean, \
            available_count[i] + 0, \
            mean(available_sum[i], available_count[i]), \
            alive_valid[i] + 0, mean(alive_sum[i], alive_valid[i]), \
            carried_count[i] + 0
    }
  }
' "$extinction_file" "$qav_file"