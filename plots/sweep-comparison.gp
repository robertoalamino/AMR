# Usage: gnuplot -c plots/sweep-comparison.gp comparison.tsv output.png
set terminal pngcairo size 1800,650
set output ARG2
set datafile separator "\t"
set multiplot layout 1,3 title "Historical and corrected sweeps: 3000 steps, 100 realisations per pair, seed 1"

unset key
set xlabel "alpha"
set ylabel "beta"
set xrange [-1.1:1.1]
set yrange [-1.1:1.1]
set xtics 0.4
set ytics 0.4
set size ratio -1

set palette defined (0 "#f2f5fa", 0.5 "#6497c5", 1 "#163b63")
set cbrange [0:1]
set cblabel "Mean final q"

set title "Historical"
plot ARG1 every ::1 using 1:2:($1-0.1):($1+0.1):($2-0.1):($2+0.1):3 \
    with boxxyerror lc palette fs solid 1.0 noborder

set title "Corrected"
plot ARG1 every ::1 using 1:2:($1-0.1):($1+0.1):($2-0.1):($2+0.1):4 \
    with boxxyerror lc palette fs solid 1.0 noborder

set palette defined (-1 "#2166ac", 0 "#ffffff", 1 "#b2182b")
set cbrange [-0.006:0.006]
set cblabel "Corrected minus historical q"
set title "Difference"
plot ARG1 every ::1 using 1:2:($1-0.1):($1+0.1):($2-0.1):($2+0.1):5 \
    with boxxyerror lc palette fs solid 1.0 noborder

unset multiplot
unset output