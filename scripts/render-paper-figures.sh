#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
if [[ $# -gt 1 ]]; then
  echo "Usage: bash scripts/render-paper-figures.sh [OUTPUT_DIRECTORY]" >&2
  exit 1
fi
if ! command -v gnuplot >/dev/null 2>&1; then
  echo "Gnuplot is required; add its installation directory to PATH." >&2
  exit 1
fi
output_directory=${1:-build/paper-figures}
mkdir -p "$output_directory/data"
# Verify source data before rendering. Never overwrite recorded result tables.
bash scripts/check-released-results.sh > "$output_directory/data-verification.log"

for implementation in historical corrected; do
  awk '
    NF!=3 {bad=1;next}
    NR>1 && $1!=previous {print ""}
    {print $0;previous=$1;rows++}
    END {if(bad || rows!=121)exit 1}
  ' "results/sweep/$implementation-aggregate.tsv" \
    > "$output_directory/data/$implementation-sweep-grid.tsv"
done

# Assert that fixed paper-style scales contain the recorded inputs.
for duration in 500 1000 2500 5000; do
  awk 'FNR==1{next} $7<0.14 || $7>0.30{bad=1} END{exit bad+0}' \
    "results/reversal/$duration/corrected-summary.tsv"
  awk 'FNR==1{next} $3<0.14 || $3>0.30{bad=1} END{exit bad+0}' \
    "results/reversal/$duration/historical-mean.tsv"
  for snapshot in last-treated final; do
    awk 'FNR==1{next} $8<0 || $8>0.4{bad=1} END{exit bad+0}' \
      "results/reversal/$duration/corrected-histogram-$snapshot.tsv"
  done
done
for implementation in historical corrected; do
  awk '$3<0 || $3>0.30{bad=1} END{exit bad+0}' \
    "results/sweep/$implementation-aggregate.tsv"
done

for format in png svg; do
  for implementation in historical corrected; do
    gnuplot -c plots/paper/sweep.gp "$implementation" \
      "$output_directory/$implementation-sweep.$format" "$format" \
      "$output_directory/data/$implementation-sweep-grid.tsv" \
      "$output_directory/data/$implementation-contours.tsv"
    gnuplot -c plots/paper/reversal-means.gp "$implementation" \
      "$output_directory/$implementation-reversal-means.$format" "$format"
  done
  for xmax in 0.6 1; do
    suffix=full
    [[ "$xmax" != 0.6 ]] || suffix=paper-range
    gnuplot -c plots/paper/reversal-distributions.gp corrected \
      "$output_directory/corrected-reversal-distributions-$suffix.$format" \
      "$format" "$xmax"
  done
done

# Check that all twelve expected figures were created in this run.
for format in png svg; do
  for figure in historical-sweep corrected-sweep \
    historical-reversal-means corrected-reversal-means \
    corrected-reversal-distributions-full corrected-reversal-distributions-paper-range; do
    [[ -s "$output_directory/$figure.$format" ]] || {
      echo "FAIL: missing or empty figure: $figure.$format" >&2
      exit 1
    }
  done
done
{
  printf 'Renderer: '
  gnuplot --version
  printf 'Data: results/checksums.sha256 (verified before rendering)\n'
  printf 'Workflow: C simulations; Bash/AWK analysis; Gnuplot figures\n'
  printf 'Distributions: fractions of all 1000 realisations; no renormalisation\n'
  printf 'Status: numerical reruns in paper-style presentation, not original published images\n'
} > "$output_directory/rendering.txt"
(
  cd "$output_directory"
  sha256sum ./*.png ./*.svg rendering.txt > figures.sha256
)
printf 'Rendered figures: %s\n' "$output_directory"
echo "Rendering completed. Inspect PNGs before committing any generated figures."
