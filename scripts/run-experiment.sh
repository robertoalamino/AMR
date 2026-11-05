#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

usage() {
  echo "Usage: bash scripts/run-experiment.sh VERSION PROGRAM SEED REALISATIONS [DURATION]"
  echo "VERSION: historical or corrected"
  echo "PROGRAM: sweep or reversal"
  echo "DURATION: required for reversal; 500, 1000, 2500 or 5000"
  exit 1
}

[[ $# -ge 4 && $# -le 5 ]] || usage

version=$1
program=$2
seed=$3
realisations=$4
duration=${5:-}

[[ "$seed" =~ ^[1-9][0-9]*$ ]] || usage
[[ "$realisations" =~ ^[1-9][0-9]*$ ]] || usage

case "$version" in
  historical)
    source_tree=historical/src
    source_label="historical comparison baseline"
    ;;
  corrected)
    source_tree=src
    source_label="corrected implementation"
    ;;
  *)
    usage
    ;;
esac

case "$program" in
  sweep)
    [[ $# -eq 4 ]] || usage
    executable=amdspace
    arguments=(--steps 3000 --realisations "$realisations" --seed "$seed")
    ;;
  reversal)
    [[ $# -eq 5 ]] || usage
    case "$duration" in
      500|1000|2500|5000) ;;
      *) usage ;;
    esac
    executable=reversaltimes
    treatment_end=$((5000 + duration))
    arguments=(
      --steps 15000
      --realisations "$realisations"
      --start 5000
      --end "$treatment_end"
      --seed "$seed"
    )
    ;;
  *)
    usage
    ;;
esac

# Source snapshots and hashes work from either a checkout or a source archive.
# Git metadata is supplementary; source hashes identify the actual run inputs.
commit=unavailable
if git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  commit=$(git rev-parse HEAD)
fi

mkdir -p build/experiments
run_dir=$(mktemp -d \
  "$PWD/build/experiments/${version}-${program}-seed${seed}.XXXXXX")

mkdir -p "$run_dir/source" "$run_dir/output"

echo "Preparing run: $run_dir"

cp -R "$source_tree" "$run_dir/source/src"
(
  cd "$run_dir/source"
  find src -type f -print0 | sort -z | xargs -0 sha256sum
) > "$run_dir/source.sha256"

gcc --version > "$run_dir/compiler.txt"
uname -a > "$run_dir/environment.txt"

sources=(
  "$run_dir/source/src/$executable.c"
  "$run_dir/source/src/genotype.c"
  "$run_dir/source/src/random.c"
  "$run_dir/source/src/response.c"
  "$run_dir/source/src/lattice.c"
  "$run_dir/source/src/config.c"
)

# This shared module exists in the corrected implementation.
if [[ -f "$run_dir/source/src/population.c" ]]; then
  sources+=("$run_dir/source/src/population.c")
fi

build_command=(
  gcc -std=c99 -Wall -Wextra -Wpedantic -O2
  "${sources[@]}"
  -o "$run_dir/$executable.exe" -lm
)

{
  printf 'Version: %s\n' "$version"
  printf 'Release checkout commit (if available): %s\n' "$commit"
  printf 'Source: %s\n' "$source_label"
  printf 'Authoritative source identification: source.sha256\n'
  printf 'Program: %s\n' "$program"
  printf 'Output directory: %s\n' "$run_dir/output"
  printf 'Build command:'
  printf ' %q' "${build_command[@]}"
  printf '\nRun command:'
  printf ' %q' "$run_dir/$executable.exe" "${arguments[@]}"
  printf '\n'
} > "$run_dir/manifest.txt"

if "${build_command[@]}" > "$run_dir/build.log" 2>&1; then
  echo "Build status: success" >> "$run_dir/manifest.txt"
else
  echo "Build status: failed" >> "$run_dir/manifest.txt"
  cat "$run_dir/build.log" >&2
  exit 1
fi

sha256sum "$run_dir/$executable.exe" \
  > "$run_dir/executable.sha256"

printf 'Started UTC: %s\n' "$(date -u +%FT%TZ)" \
  >> "$run_dir/manifest.txt"

echo "Running $version $program..."
TIMEFORMAT='Elapsed seconds: %3R'

if {
  time (
    cd "$run_dir/output"
    "$run_dir/$executable.exe" "${arguments[@]}"
  )
} > "$run_dir/stdout.log" 2> "$run_dir/run.log"
then
  status=0
else
  status=$?
fi

{
  printf 'Finished UTC: %s\n' "$(date -u +%FT%TZ)"
  printf 'Exit status: %s\n' "$status"
} >> "$run_dir/manifest.txt"

(
  cd "$run_dir/output"
  find . -type f -exec sha256sum {} +
) > "$run_dir/output.sha256"

cat "$run_dir/run.log"
echo "Run records: $run_dir"

if [[ "$status" -ne 0 ]]; then
  echo "FAIL: simulation exited with status $status." >&2
  exit "$status"
fi

echo "Run completed; scientific output validation is still required."
