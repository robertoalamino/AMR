#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."
repository_root="$PWD"

test_dir=$(mktemp -d "$PWD/build/extinction-during.XXXXXX")
trap 'rm -rf -- "$test_dir"' EXIT

# Create a temporary test copy that forces extinction just before
# the culling stage of step 3. The production source is untouched.
awk '
  /\/\/ Culling cell population if necessary/ {
    matches++
    print "      if(t == 3) {"
    print "        if(no == 0) {"
    print "          fprintf(stderr, \"Test population extinct too early\\n\");"
    print "          return EXIT_FAILURE;"
    print "        }"
    print "        for(int test_i = 1; test_i <= no; test_i++)"
    print "          s[oc[test_i]] = 0;"
    print "        nno = 0;"
    print "      }"
  }
  { print }
  END {
    if(matches != 1) {
      print "Expected one culling insertion point; found " matches > "/dev/stderr"
      exit 1
    }
  }
' src/reversaltimes.c > "$test_dir/reversaltimes-test.c"

gcc -std=c99 -Wall -Wextra -Wpedantic -O2 -Isrc \
  "$test_dir/reversaltimes-test.c" \
  src/genotype.c src/random.c src/response.c \
  src/lattice.c src/config.c src/population.c \
  -o "$test_dir/reversaltimes-test.exe" -lm

(
  cd "$test_dir"
  ./reversaltimes-test.exe \
    --steps 6 --realisations 1 --start 2 --end 4 --seed 1
)

# The first extinction is the deliberately forced event at step 3.
awk '
  NR != 1 || NF != 2 || $1 != 1 || $2 != 3 { bad = 1 }
  END {
    if(bad || NR != 1) {
      print "FAIL: expected extinction at step 3."
      exit 1
    }
  }
' "$test_dir/extinction.txt"

# qav columns represent steps 2 through 6. The value measured at
# step 3 must remain unchanged at steps 4, 5 and 6.
awk '
  NR != 1 || NF != 5 { bad = 1 }
  {
    if($2 !~ /^[0-9]+([.][0-9]+)?$/)
      bad = 1
    if($3 != $2 || $4 != $2 || $5 != $2)
      bad = 1
  }
  END {
    if(bad || NR != 1) {
      print "FAIL: last measured susceptibility was not carried forward."
      exit 1
    }
  }
' "$test_dir/qav.txt"

echo "PASS: extinction time and carried-forward susceptibility."