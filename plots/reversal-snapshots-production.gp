# Usage:
# gnuplot -c plots/reversal-snapshots-production.gp FIRST.tsv SECOND.tsv OUTPUT.png DURATION

set terminal pngcairo size 1400,900 font "Arial,12"
set output ARG3

set datafile separator "\t"
set datafile missing "NaN"

set multiplot layout 2,1 title \
    sprintf("Corrected reversal — %s-step treatment, 1000 realisations, seed 1", ARG4)

set lmargin 12
set rmargin 3
set xrange [0:1]
set yrange [0:*]
set xlabel "Recorded death probability"
set ylabel "Fraction of all realisations"
set boxwidth 0.01 absolute
set style fill solid 1.0 noborder
set grid ytics
set key outside top center horizontal

# Column 8 is the available bin count divided by all realisations.
# Draw cumulative heights, largest first, to show stacked categories.
set title "Last treated step"
plot ARG1 every ::1 using 3:8 \
         with boxes lc rgb "#0072B2" title "Alive after step", \
     ARG1 every ::1 using 3:($4 > 0 ? ($6+$7)*$8/$4 : 0) \
         with boxes lc rgb "#D55E00" title "Extinct before step", \
     ARG1 every ::1 using 3:($4 > 0 ? $7*$8/$4 : 0) \
         with boxes lc rgb "#009E73" title "Extinction at step"

set title "Final recorded step"
plot ARG2 every ::1 using 3:8 \
         with boxes lc rgb "#0072B2" title "Alive after step", \
     ARG2 every ::1 using 3:($4 > 0 ? ($6+$7)*$8/$4 : 0) \
         with boxes lc rgb "#D55E00" title "Extinct before step", \
     ARG2 every ::1 using 3:($4 > 0 ? $7*$8/$4 : 0) \
         with boxes lc rgb "#009E73" title "Extinction at step"

unset multiplot