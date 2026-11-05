# Preserve -c arguments before loading shared styles.
implementation=ARG1
output_file=ARG2
output_format=ARG3
display_limit=ARG4

# gnuplot -c plots/paper/reversal-distributions.gp corrected OUTPUT FORMAT XMAX
# No historical individual snapshots are included in the release data.
if (implementation ne "corrected") exit error "Only corrected distribution data are available"
xmax=display_limit+0
if (!(xmax==0.6 || xmax==1)) exit error "Expected xmax 0.6 or 1"
image_width=1500
image_height=1100
load "plots/paper/style.gp"

set multiplot layout 2,2 rowsfirst title "Corrected reversal distributions: paper-style view"
set xrange [0:xmax]
# A common scale covers the largest released bin, without clipping.
set yrange [0:0.40]
set xtics 0.1
set ytics 0.1
set format x "%.1f"
set format y "%.2f"
set lmargin 9
set rmargin 3
set bmargin 4
set tmargin 3
set key top left horizontal font ",11" spacing 1.2
set style fill transparent solid 0.65 noborder

array durations[4]=[5000,2500,1000,500]
do for [panel=1:4] {
    duration=durations[panel]
    last=sprintf("results/reversal/%d/corrected-histogram-last-treated.tsv",duration)
    final=sprintf("results/reversal/%d/corrected-histogram-final.tsv",duration)
    # Stats must not inherit the display limits and filter tail counts.
    set xrange [*:*]
    set yrange [*:*]
    stats last every ::1 using ($1>=xmax ? $4 : 0) nooutput
    last_outside=STATS_sum
    stats final every ::1 using ($1>=xmax ? $4 : 0) nooutput
    final_outside=STATS_sum
    set xrange [0:xmax]
    set yrange [0:0.40]
    # q = 1 is inside the final bin in the full-range view.
    set label 1 sprintf("T* = %d",duration) at graph 0.96,0.87 right
    if (xmax<1) {
        set label 2 sprintf("q >= 0.6: last %d/1000; final %d/1000",int(last_outside),int(final_outside)) at graph 0.96,0.78 right font ",10"
    } else {
        unset label 2
    }
    set title ""
    plot last every ::1 using ($3<xmax ? $3 : 1/0):8 with filledcurves y1=0 lc rgb "#4472c4" title "Last treated step", \
         final every ::1 using ($3<xmax ? $3 : 1/0):8 with filledcurves y1=0 lc rgb "#ffe6a3" title "Final step", \
         last every ::1 using ($3<xmax ? $3 : 1/0):8 with lines lw 0.8 lc rgb "#1f355a" notitle, \
         final every ::1 using ($3<xmax ? $3 : 1/0):8 with lines lw 0.8 lc rgb "#756345" notitle
}
unset multiplot
unset output
