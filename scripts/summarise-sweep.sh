#!/usr/bin/env bash
set -euo pipefail

# Read sweep-realisations.txt from standard input; write TSV to stdout.
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

  function mean(sum, count) {
    return count ? sprintf("%.17g", sum / count) : "NaN"
  }

  {
    if (NF != 5 || !is_number($1) || !is_number($2) ||
        $3 !~ /^[0-9]+$/ || $3 < 1 ||
        $5 !~ /^-?[0-9]+$/ || $5 < -1) {
      print "Invalid record on line " NR > "/dev/stderr"
      bad = 1
      next
    }

    alpha = $1
    beta = $2
    value = $4
    extinction_time = $5 + 0
    key = alpha SUBSEP beta

    if (!(key in seen)) {
      seen[key] = 1
      order[++pair_count] = key
      pair_alpha[key] = alpha
      pair_beta[key] = beta
    }

    total[key]++

    if (extinction_time == -1)
      alive[key]++
    else
      extinct[key]++

    if (is_nan(value)) {
      unavailable[key]++
    } else if (!is_number(value) || value < 0 || value > 1) {
      print "Invalid susceptibility on line " NR > "/dev/stderr"
      bad = 1
    } else {
      available_count[key]++
      available_sum[key] += value

      if (extinction_time == -1) {
        alive_valid[key]++
        alive_sum[key] += value
      }
    }
  }

  END {
    if (NR == 0) {
      print "No realisation records supplied." > "/dev/stderr"
      exit 1
    }

    if (bad)
      exit 1

    print "alpha", "beta", "total", "alive", "extinct", \
          "unavailable", "historical_style_mean", \
          "available_count", "available_mean", \
          "alive_valid_count", "alive_mean"

    for (i = 1; i <= pair_count; i++) {
      key = order[i]

      # The C aggregate propagates NaN if any realisation is unavailable.
      historical_mean = unavailable[key] \
        ? "NaN" : mean(available_sum[key], available_count[key])

      print pair_alpha[key], pair_beta[key], total[key], \
            alive[key] + 0, extinct[key] + 0, unavailable[key] + 0, \
            historical_mean, available_count[key] + 0, \
            mean(available_sum[key], available_count[key]), \
            alive_valid[key] + 0, mean(alive_sum[key], alive_valid[key])
    }
  }
'