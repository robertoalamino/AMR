# Preserve -c arguments before loading shared styles.
implementation=ARG1
output_file=ARG2
output_format=ARG3

# gnuplot -c plots/paper/reversal-means.gp IMPLEMENTATION OUTPUT FORMAT
if (!(implementation eq "corrected" || implementation eq "historical")) exit error "Invalid implementation"
image_width=1200
image_height=800
load "plots/paper/style.gp"

set title sprintf("%s reversal: paper-style view", implementation eq "corrected" ? "Corrected" : "Historical baseline")
set xlabel "Steps since treatment onset"
set ylabel "Mean recorded q"
set xrange [0:10000]
# These limits cover every mean in the verified released inputs.
set yrange [0.14:0.30]
set xtics 2000
set ytics 0.02
set lmargin 10
set rmargin 3
set bmargin 5
set key bottom right title "Treatment duration (steps)" opaque font ",12"
set label 1 "1000 realisations; m = 0.001; seed 1; means include carried observations" at screen 0.5,0.025 center font ",11"
file(d)=sprintf("results/reversal/%d/%s", d, implementation eq "corrected" ? "corrected-summary.tsv" : "historical-mean.tsv")
column_q=(implementation eq "corrected" ? 7 : 3)
plot file(500) every ::1 using ($1-5000):(column(column_q)) with lines lw 2 lc rgb "#e5aa00" title "500", \
     file(1000) every ::1 using ($1-5000):(column(column_q)) with lines lw 2 lc rgb "#528535" title "1000", \
     file(2500) every ::1 using ($1-5000):(column(column_q)) with lines lw 2 lc rgb "#dd2626" title "2500", \
     file(5000) every ::1 using ($1-5000):(column(column_q)) with lines lw 2 lc rgb "#111111" title "5000"
unset output
