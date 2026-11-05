#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

test_dir=$(mktemp -d "$PWD/build/amdspace-extinction-during.XXXXXX")
trap 'rm -rf -- "$test_dir"' EXIT

# Force extinction at step 3, after the death-probability measurement.
# Record that measurement independently for comparison with the final output.
awk '
  /\/\/ Culling cell population if necessary/ {
    matches++
    print "          if(t == 3) {"
    print "            if(no == 0 || !isfinite(apt[t])) {"
    print "              fprintf(stderr, \"No valid measurement at step 3\\n\");"
    print "              return EXIT_FAILURE;"
    print "            }"
    print "            FILE *test_expected = fopen(\"expected.txt\", \"a\");"
    print "            if(test_expected == NULL) {"
    print "              perror(\"Cannot open expected.txt\");"
    print "              return EXIT_FAILURE;"
    print "            }"
    print "            fprintf(test_expected, \"%f\\t%f\\t%f\\n\", aa, ab, apt[t]);"
    print "            if(fclose(test_expected) == EOF) {"
    print "              perror(\"Cannot finish expected.txt\");"
    print "              return EXIT_FAILURE;"
    print "            }"
    print "            for(int test_i = 1; test_i <= no; test_i++)"
    print "              s[oc[test_i]] = 0;"
    print "            nno = 0;"
    print "          }"
  }
  { print }
  END {
    if(matches != 1) {
      print "Expected one culling insertion point; found " matches > "/dev/stderr"
      exit 1
    }
  }
' src/amdspace.c > "$test_dir/amdspace-test.c"

gcc -std=c99 -Wall -Wextra -Wpedantic -O2 -Isrc \
  "$test_dir/amdspace-test.c" \
  src/genotype.c src/random.c src/response.c \
  src/lattice.c src/config.c src/population.c \
  -o "$test_dir/amdspace-test.exe" -lm

(
  cd "$test_dir"
  ./amdspace-test.exe --steps 5 --realisations 1 --seed 1
)

# Each of the 121 parameter combinations has one realisation,
# deliberately made extinct at step 3.
awk '
  NF != 4 || $3 != 1 || $4 != 3 { bad = 1 }
  END {
    if(bad || NR != 121) {
      print "FAIL: expected extinction at step 3 for every combination."
      exit 1
    }
  }
' "$test_dir/extinction.txt"

# With one realisation, the final reported value must equal the
# measurement taken immediately before forced extinction.
cmp "$test_dir/expected.txt" "$test_dir/asymdrug.txt"

echo "PASS: sweep extinction time and carried-forward measurement."