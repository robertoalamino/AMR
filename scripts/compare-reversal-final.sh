#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

if [[ $# -ne 2 ]]; then
  echo "Usage: bash scripts/compare-reversal-final.sh HISTORICAL_RUN CORRECTED_RUN" >&2
  exit 1
fi

awk '
  function fail(message) {
    print "FAIL: " message > "/dev/stderr"
    bad = 1
  }

  # Stable accumulation of means and variances.
  function add_value(key, value, delta) {
    n[key]++
    delta = value - mean[key]
    mean[key] += delta/n[key]
    m2[key] += delta*(value-mean[key])
  }

  function report_mean(label, key, h, c, se, difference) {
    h = "1:" key
    c = "2:" key
    if (n[h] < 2 || n[c] < 2) {
      fail("insufficient observations for " label)
      return
    }

    se = sqrt(m2[h]/(n[h]-1)/n[h] + m2[c]/(n[c]-1)/n[c])
    difference = mean[c]-mean[h]

    printf "%s\n", label
    printf "  Historical: %.9f (n=%d)\n", mean[h], n[h]
    printf "  Corrected:  %.9f (n=%d)\n", mean[c], n[c]
    printf "  Corrected minus historical: %.9f\n", difference
    printf "  Approximate 95%% interval: [%.9f, %.9f]\n\n",
           difference-1.96*se, difference+1.96*se
  }

  # Wilson intervals for individual proportions.
  function wilson(successes, total, key, p, denominator, centre, half) {
    p = successes/total
    denominator = 1+1.96^2/total
    centre = (p+1.96^2/(2*total))/denominator
    half = 1.96*sqrt(p*(1-p)/total + 1.96^2/(4*total^2))/denominator
    lower[key] = centre-half
    upper[key] = centre+half
  }

  # Newcombe interval for a difference of independent proportions.
  function report_fraction(label, hs, hn, cs, cn, ph, pc, d, lo, hi) {
    if (hn == 0 || cn == 0) {
      fail("empty denominator for " label)
      return
    }

    ph = hs/hn
    pc = cs/cn
    wilson(hs, hn, "h")
    wilson(cs, cn, "c")
    d = pc-ph

    lo = d-sqrt((pc-lower["c"])^2+(upper["h"]-ph)^2)
    hi = d+sqrt((upper["c"]-pc)^2+(ph-lower["h"])^2)

    printf "%s\n", label
    printf "  Historical: %d/%d = %.2f%%\n", hs, hn, 100*ph
    printf "  Corrected:  %d/%d = %.2f%%\n", cs, cn, 100*pc
    printf "  Difference: %.2f percentage points\n", 100*d
    printf "  Approximate 95%% interval: [%.2f, %.2f] percentage points\n\n",
           100*lo, 100*hi
  }

  FILENAME == ARGV[1] || FILENAME == ARGV[3] {
    group = (FILENAME == ARGV[1]) ? 1 : 2
    if (NF != 15001)
      fail("unexpected population dimensions")

    alive[group,FNR] = ($NF > 0)
    population_rows[group]++
    if ($NF == 0) extinct[group]++
    next
  }

  {
    group = (FILENAME == ARGV[2]) ? 1 : 2
    q_rows[group]++

    if (NF != 10001) {
      fail("unexpected susceptibility dimensions")
      next
    }

    value = $NF
    if (value !~ /^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/ ||
        value < 0 || value > 1) {
      fail("unavailable or invalid final susceptibility")
      next
    }

    add_value(group ":all", value)
    if (value <= 0.1) low_all[group]++
    if (value >= 0.27) high_all[group]++

    if (alive[group,FNR]) {
      add_value(group ":alive", value)
      if (value <= 0.1) low_alive[group]++
      if (value >= 0.27) high_alive[group]++
    }
  }

  END {
    for (group = 1; group <= 2; group++)
      if (population_rows[group] != 1000 || q_rows[group] != 1000)
        fail("expected 1000 realisations per implementation")

    if (bad) exit 1

    report_mean("Final q including carried values", "all")
    report_mean("Final q among final recorded survivors", "alive")

    report_fraction("Zero final recorded population",
      extinct[1], population_rows[1],
      extinct[2], population_rows[2])

    report_fraction("Final q >= 0.27, all realisations",
      high_all[1], n["1:all"], high_all[2], n["2:all"])

    report_fraction("Final q <= 0.1, all realisations",
      low_all[1], n["1:all"], low_all[2], n["2:all"])

    report_fraction("Final q >= 0.27, survivors",
      high_alive[1], n["1:alive"], high_alive[2], n["2:alive"])

    report_fraction("Final q <= 0.1, survivors",
      low_alive[1], n["1:alive"], low_alive[2], n["2:alive"])

    if (bad) exit 1
  }
' "$1/output/population.txt" "$1/output/qav.txt" \
  "$2/output/population.txt" "$2/output/qav.txt"