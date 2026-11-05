# Usage: gnuplot -c plots/sweep-outcomes.gp input.tsv output.png "Title"
set terminal pngcairo size 1400,650
set output ARG2
set datafile separator "\t"
set datafile missing "NaN"

set multiplot layout 1,2 title ARG3
unset key
set xlabel "alpha"
set ylabel "beta"
set xrange [-1.1:1.1]
set yrange [-1.1:1.1]
set xtics 0.4
set ytics 0.4
set size ratio -1
set cbrange [0:1]
set palette defined (0 "#f2f5fa", 0.5 "#6497c5", 1 "#163b63")

set title "Mean final susceptibility among survivors"
set cblabel "Mean final q"
plot ARG1 every ::1 using 1:2:($1-0.1):($1+0.1):($2-0.1):($2+0.1):3 \
    with boxxyerror lc palette fs solid 1.0 noborder

set title "Fraction extinct by end of run"
set cblabel "Extinction fraction"
plot ARG1 every ::1 using 1:2:($1-0.1):($1+0.1):($2-0.1):($2+0.1):4 \
    with boxxyerror lc palette fs solid 1.0 noborder

unset multiplot
unset output