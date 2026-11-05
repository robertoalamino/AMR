#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

mkdir -p build
test_dir=$(mktemp -d "$PWD/build/parent-selection.XXXXXX")
trap 'rm -rf -- "$test_dir"' EXIT

for program in reversaltimes amdspace
do
  awk '
    /\/\/ CELL REPRODUCTION/ {
      reproduction = 1
      starts++
      print
      print "int test_seen[N+1];"
      print "for(int k = 0; k <= N; k++) test_seen[k] = 0;"
      print "r = atof(getenv(\"AMR_TEST_RATE\"));"
      next
    }

    reproduction && /iaux = .*ran2/ {
      # Always select the first remaining parent. This makes a skipped
      # removal produce a deterministic duplicate on the next attempt.
      choices++
      print "iaux = 1;"
      next
    }

    reproduction && /cc = oc\[iaux\];/ {
      selections++
      print
      print "if(cc < 1 || cc > N || ++test_seen[cc] != 1) {"
      print "  fprintf(stderr, \"FAIL: repeated or invalid parent\\n\");"
      print "  return EXIT_FAILURE;"
      print "}"
      next
    }

    reproduction && /no = nno;.*Updating/ {
      finishes++
      print "if(no == 0) {"
      print "  fprintf(stderr, \"FAIL: test has no parents\\n\");"
      print "  return EXIT_FAILURE;"
      print "}"
      print "for(int k = 1; k <= no; k++) {"
      print "  if(test_seen[oc[k]] != 1) {"
      print "    fprintf(stderr, \"FAIL: original parent was missed\\n\");"
      print "    return EXIT_FAILURE;"
      print "  }"
      print "}"
      print "fprintf(stderr, \"PASS: each original parent considered once\\n\");"
      reproduction = 0
    }

    { print }

    END {
      if(starts != 1 || choices != 1 ||
         selections != 1 || finishes != 1) {
        print "FAIL: unexpected reproduction-loop structure" > "/dev/stderr"
        exit 1
      }
    }
  ' "src/$program.c" > "$test_dir/$program-test.c"

  gcc -std=c99 -Wall -Wextra -Wpedantic -O2 -Isrc \
    "$test_dir/$program-test.c" \
    src/genotype.c src/random.c src/response.c \
    src/lattice.c src/config.c src/population.c \
    -o "$test_dir/$program-test.exe" -lm

  for rate in 0 0.5 1
  do
    run_dir="$test_dir/$program-$rate"
    mkdir -p "$run_dir"

    (
      cd "$run_dir"

      if [[ "$program" == reversaltimes ]]; then
        AMR_TEST_RATE="$rate" ../reversaltimes-test.exe \
          --steps 1 --realisations 1 --start 1 --end 2 --seed 1
      else
        AMR_TEST_RATE="$rate" ../amdspace-test.exe \
          --steps 1 --realisations 1 --seed 1
      fi
    ) > "$run_dir/stdout.txt" 2> "$run_dir/stderr.txt" || {
      cat "$run_dir/stderr.txt" >&2
      exit 1
    }

    expected=1
    if [[ "$program" == amdspace ]]; then
      expected=121
    fi

    count=$(grep -c '^PASS: each original parent considered once$' \
      "$run_dir/stderr.txt" || true)

    if [[ "$count" != "$expected" ]]; then
      cat "$run_dir/stderr.txt" >&2
      echo "FAIL: unexpected number of parent-selection checks." >&2
      exit 1
    fi

    echo "PASS: $program parent selection with r=$rate."
  done
done