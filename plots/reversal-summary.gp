# Arguments: summary.tsv output.png treatment_start treatment_end

input_file = ARG1
output_file = ARG2
treatment_start = ARG3 + 0
treatment_end = ARG4 + 0
duration = treatment_end - treatment_start

set terminal pngcairo size 1400,1000 font "Arial,12"
set output output_file

set datafile separator "\t"
set datafile columnheaders
set datafile missing "NaN"

set border 3
set tics nomirror
set grid back lc rgb "#dddddd"
set key outside top center horizontal
set lmargin 12
set rmargin 3
set xrange [0:*]

set multiplot layout 2,1 title \
    "Short reversal demonstration — plotting check, not published-result reproduction"

set xlabel "Steps since treatment onset"
set ylabel "Recorded death probability"
set yrange [*:*]

set arrow 1 from first duration, graph 0 \
             to first duration, graph 1 \
             nohead dt 2 lc rgb "#777777"

plot input_file using ($1-treatment_start):7 \
       with lines lw 2 lc rgb "#333333" \
       title "All: NaN propagates", \
     input_file using ($1-treatment_start):9 \
       with lines lw 2 dt 2 lc rgb "#d97706" \
       title "Available: includes carry-forward", \
     input_file using ($1-treatment_start):11 \
       with lines lw 2 dt 3 lc rgb "#0072b2" \
       title "Recorded values: alive after step"

set ylabel "Fraction of realisations"
set yrange [-0.02:1.05]

plot input_file using ($1-treatment_start):($3/$2) \
       with lines lw 2 lc rgb "#0072b2" \
       title "Alive at end of step", \
     input_file using ($1-treatment_start):(($4+$5)/$2) \
       with lines lw 2 lc rgb "#cc3311" \
       title "Extinct by end of step", \
     input_file using ($1-treatment_start):($8/$2) \
       with lines lw 2 dt 2 lc rgb "#228833" \
       title "Available susceptibility"

unset multiplot
unset output