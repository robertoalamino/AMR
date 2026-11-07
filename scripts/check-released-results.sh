#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
(cd results && sha256sum --strict -c checksums.sha256)
mkdir -p build
check_dir=$(mktemp -d "$PWD/build/released-results.XXXXXX")
trap 'rm -rf -- "$check_dir"' EXIT

# Require identical table shape/labels and compare numeric fields with a
# tolerance below the six-decimal simulation output precision.
compare_tables() {
  awk '
    function number(x) { return x ~ /^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/ }
    FILENAME == ARGV[1] {
      rows=FNR; width[FNR]=NF
      for(i=1;i<=NF;i++) value[FNR,i]=$i
      next
    }
    {
      actual_rows=FNR
      if(FNR>rows || NF!=width[FNR]) bad=1
      for(i=1;i<=NF;i++) {
        a=value[FNR,i]; b=$i
        if(number(a) && number(b)) {
          d=a-b; if(d<0) d=-d
          if(d>1e-12) bad=1
        } else if(a!=b) bad=1
      }
    }
    END { if(bad || rows!=actual_rows || !rows) exit 1 }
  ' "$1" "$2" || {
    echo "FAIL: regenerated table differs: $1" >&2
    return 1
  }
}

for duration in 500 1000 2500 5000; do
  directory="results/reversal/$duration"
  awk '
    FNR==1 { next }
    {
      if(NF!=3 || $1!=4998+FNR || $2!=FNR-2 ||
         $3 !~ /^[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)?$/ || $3<0 || $3>1) bad=1
      rows++
    }
    END { if(bad || rows!=10001) exit 1 }
  ' "$directory/historical-mean.tsv"
  for kind in last-treated withdrawal final; do
    snapshot="$directory/corrected-snapshot-$kind.tsv"
    step=15000
    [[ "$kind" != last-treated ]] || step=$((4999+duration))
    [[ "$kind" != withdrawal ]] || step=$((5000+duration))
    awk -v step="$step" '
      FILENAME == ARGV[1] {
        if(NF!=2 || $1!=FNR || $2 !~ /^(-1|[0-9]+)$/ || $2>15000) bad=1
        time[$1]=$2; n++; next
      }
      FNR==1 { next }
      {
        rows++; id=$1
        if(NF!=6 || id!=rows || $2!=step || !(id in time) || $4!=time[id]) bad=1
        expected=(time[id]==-1 || time[id]>step) ? "alive" : (time[id]==step ? "extinction_step" : "extinct_before")
        if($5!=expected) bad=1
      }
      END { if(bad || n!=1000 || rows!=1000) exit 1 }
    ' "$directory/corrected-extinction.tsv" "$snapshot"
    bash scripts/bin-reversal-snapshot.sh < "$snapshot" \
      > "$check_dir/histogram.tsv" 2> "$check_dir/counts.txt"
    compare_tables "$directory/corrected-histogram-$kind.tsv" "$check_dir/histogram.tsv"
    awk -v step="$step" '
      FILENAME==ARGV[1] {
        if(FNR==1) next
        if($1!=4998+FNR || $2!=1000 || $3+$4+$5!=1000) bad=1
        if($1==step) { expected=$7; found++ }
        summary_rows++; next
      }
      FNR==1 { next }
      { if($6==1) { sum+=$3; n++ } }
      END {
        difference=sum/n-expected; if(difference<0) difference=-difference
        if(bad || summary_rows!=10001 || found!=1 || n!=1000 || difference>1e-12) exit 1
      }
    ' "$directory/corrected-summary.tsv" "$snapshot"
  done
  echo "PASS: $duration-step snapshots, extinction associations, means and regenerated histograms."
done

bash scripts/summarise-sweep.sh < results/sweep/corrected-realisations.tsv \
  > "$check_dir/sweep-summary.tsv"
compare_tables results/sweep/corrected-summary.tsv "$check_dir/sweep-summary.tsv"
awk '
  FILENAME==ARGV[1] {
    key=$1 SUBSEP $2 SUBSEP $3
    if(NF!=4 || key in extinction) bad=1
    extinction[key]=$4; n++; next
  }
  {
    key=$1 SUBSEP $2 SUBSEP $3
    if(NF!=5 || !(key in extinction) || key in seen || $5!=extinction[key]) bad=1
    seen[key]=1; rows++
  }
  END { if(bad || n!=12100 || rows!=12100) exit 1 }
' results/sweep/corrected-extinction.tsv results/sweep/corrected-realisations.tsv
awk '
  FILENAME==ARGV[1] { a[$1,$2]=$3; count++; next }
  FNR==1 {next}
  {
    key=$1 SUBSEP $2
    if(!(key in a) || key in seen || NF!=11 || $3!=100 || $4+$5!=$3) bad=1
    seen[key]=1; rows++; d=$7-a[key]; if(d<0)d=-d
    if(d>0.00000051) bad=1
  }
  END { if(bad || count!=121 || rows!=121) exit 1 }
' results/sweep/corrected-aggregate.tsv results/sweep/corrected-summary.tsv
for version in historical corrected; do
  column=3
  [[ "$version" != corrected ]] || column=4
  awk -v column="$column" '
    FILENAME==ARGV[1] {a[$1,$2]=$3; count++;next}
    FNR==1 {next}
    {
      key=$1 SUBSEP $2
      if(!(key in a) || key in seen || NF!=5)bad=1
      seen[key]=1;rows++;d=$column-a[key];if(d<0)d=-d
      if(d>1e-12)bad=1
      d=($4-$3)-$5;if(d<0)d=-d;if(d>1e-12)bad=1
    }
    END {if(bad || count!=121 || rows!=121)exit 1}
  ' "results/sweep/$version-aggregate.tsv" results/sweep/comparison.tsv
done
echo "PASS: sweep individual records, regenerated summary, aggregates and comparison."
echo "PASS: all released-result checks completed."
