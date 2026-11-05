#!/usr/bin/env bash
set -euo pipefail

# Read an extracted snapshot from standard input; write histogram TSV.
awk '
  BEGIN {
    OFS = "\t"
  }

  function is_nan(value) {
    return tolower(value) ~ /^[-+]?nan(\([^)]*\))?$/
  }

  function is_number(value) {
    return value ~ /^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/
  }

  function fraction(count, denominator) {
    return denominator ? sprintf("%.17g", count / denominator) : "NaN"
  }

  function fail(message) {
    print "FAIL: " message > "/dev/stderr"
    bad = 1
  }

  NR == 1 {
    if (NF != 6 || $1 != "realisation" || $3 != "susceptibility")
      fail("unexpected snapshot header")
    next
  }

  {
    if (NF != 6 ||
        ($5 != "alive" && $5 != "extinct_before" &&
         $5 != "extinction_step")) {
      fail("invalid snapshot row " NR)
      next
    }

    if (total == 0)
      snapshot_step = $2
    else if ($2 != snapshot_step)
      fail("snapshot contains multiple observation steps")

    total++

    if (is_nan($3)) {
      if ($6 != 0)
        fail("NaN marked as available on row " NR)
      missing++
      next
    }

    if ($6 != 1 || !is_number($3) || $3 < 0 || $3 > 1) {
      fail("invalid available measurement on row " NR)
      next
    }

    value = $3 + 0

    # Compare against edges directly, avoiding multiplication followed
    # by truncation near decimal boundaries. Include 1 in the last bin.
    bin = 0
    while (bin < 99 && value >= (bin + 1) / 100)
      bin++

    available_bins[bin]++
    available_total++

    if ($5 == "alive") {
      alive_bins[bin]++
      alive_available++
    } else if ($5 == "extinct_before") {
      extinct_bins[bin]++
    } else {
      extinction_step_bins[bin]++
    }
  }

  END {
    if (total == 0)
      fail("snapshot contains no realisations")

    if (bad)
      exit 1

    print "step=" snapshot_step, "total=" total, \
          "available=" available_total + 0, \
          "unavailable=" missing + 0, \
          "alive_available=" alive_available + 0 > "/dev/stderr"

    print "lower", "upper", "centre", "available_count", \
          "alive_count", "extinct_before_count", \
          "extinction_step_count", "fraction_total", \
          "fraction_available", "fraction_alive_available"

    for (bin = 0; bin < 100; bin++) {
      print sprintf("%.2f", bin / 100), \
            sprintf("%.2f", (bin + 1) / 100), \
            sprintf("%.3f", (bin + 0.5) / 100), \
            available_bins[bin] + 0, alive_bins[bin] + 0, \
            extinct_bins[bin] + 0, extinction_step_bins[bin] + 0, \
            fraction(available_bins[bin], total), \
            fraction(available_bins[bin], available_total), \
            fraction(alive_bins[bin], alive_available)
    }
  }
'