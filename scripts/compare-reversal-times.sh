#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

if [[ $# -ne 4 ]]; then
  echo "Usage: bash scripts/compare-reversal-times.sh HISTORICAL_RUN CORRECTED_RUN ANALYSIS_DIR DURATION" >&2
  exit 1
fi

historical_dir=$1
corrected_dir=$2
analysis_dir=$3
duration=$4

case "$duration" in
  500|1000|2500|5000) ;;
  *)
    echo "Duration must be 500, 1000, 2500 or 5000." >&2
    exit 1
    ;;
esac

mkdir -p "$analysis_dir"

# Settings shared by the documented reversal comparisons.
start=5000
end=$((start + duration))
steps=15000
realisations=1000

for version in historical corrected
do
  if [[ "$version" == historical ]]; then
    input_dir=$historical_dir
  else
    input_dir=$corrected_dir
  fi

  awk -v start="$start" -v steps="$steps" \
      -v expected_rows="$realisations" '
    BEGIN {
      OFS = "\t"
      expected_columns = steps-start+1
    }

    {
      if (NF != expected_columns) {
        print "FAIL: unexpected susceptibility dimensions on row " NR \
          > "/dev/stderr"
        bad = 1
        next
      }

      for (i = 1; i <= NF; i++) {
        value = $i

        if (value !~ /^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/ ||
            value < 0 || value > 1) {
          print "FAIL: invalid or unavailable susceptibility on row " NR \
            > "/dev/stderr"
          bad = 1
          exit
        }

        sum[i] += value
      }
    }

    END {
      if (NR != expected_rows) {
        print "FAIL: unexpected realisation count" > "/dev/stderr"
        bad = 1
      }
      if (bad) exit 1

      print "step", "elapsed", "mean_q"

      for (i = 1; i <= expected_columns; i++)
        printf "%d\t%d\t%.17g\n",
               start+i-1, i-1, sum[i]/NR
    }
  ' "$input_dir/output/qav.txt" \
    > "$analysis_dir/$version-mean.tsv"
done

awk -v start="$start" -v end="$end" -v steps="$steps" '
  function analyse(group, label, j, fraction, amplitude,
                   acquisition_target, recovery_target,
                   acquisition, recovery, t) {
    amplitude = q[group,start]-q[group,end-1]

    printf "\n%s: means including carried values\n", label
    printf "Onset mean: %.9f; last-treated mean: %.9f\n",
           q[group,start], q[group,end-1]

    if (amplitude <= 0) {
      print "FAIL: no positive treatment-associated decrease."
      bad = 1
      return
    }

    for (j = 1; j <= 2; j++) {
      fraction = (j == 1) ? 0.5 : 0.9
      acquisition_target = q[group,start]-fraction*amplitude
      recovery_target = q[group,end-1]+fraction*amplitude

      acquisition = -1
      recovery = -1

      for (t = start; t < end; t++)
        if (q[group,t] <= acquisition_target) {
          acquisition = t-start
          break
        }

      for (t = end; t <= steps; t++)
        if (q[group,t] >= recovery_target) {
          recovery = t-end
          break
        }

      printf "%.0f%% of treatment-associated change:\n", 100*fraction

      if (acquisition < 0)
        print "  Acquisition crossing: not reached"
      else
        printf "  Acquisition crossing: %d steps since onset\n", acquisition

      if (recovery < 0)
        print "  Recovery crossing: not reached within observation window"
      else
        printf "  Recovery crossing: %d steps since withdrawal\n", recovery
    }
  }

  FNR == 1 { next }

  {
    group = (FILENAME == ARGV[1]) ? 1 : 2
    q[group,$1] = $3
    rows[group]++
  }

  END {
    expected = steps-start+1
    if (rows[1] != expected || rows[2] != expected) {
      print "FAIL: incomplete mean trajectories." > "/dev/stderr"
      exit 1
    }

    print "Treatment interval: [" start ", " end ")"
    print "Final step: " steps
    print "Each implementation uses its own treatment-associated amplitude."
    print "Times are first crossings; crossing-time uncertainty is not estimated."

    analyse(1, "Historical")
    analyse(2, "Corrected")

    if (bad) exit 1
  }
' "$analysis_dir/historical-mean.tsv" \
  "$analysis_dir/corrected-mean.tsv" |
  tee "$analysis_dir/crossing-times.txt"

{
  printf 'Analysis commit: %s\n' "$(git rev-parse HEAD)"
  printf 'Historical run: %s\n' "$historical_dir"
  printf 'Corrected run: %s\n' "$corrected_dir"
  printf 'Settings: start=%s end=%s steps=%s realisations=%s\n' \
    "$start" "$end" "$steps" "$realisations"
  sha256sum scripts/compare-reversal-times.sh \
    "$historical_dir/output/qav.txt" \
    "$corrected_dir/output/qav.txt"
} > "$analysis_dir/crossing-times-provenance.txt"