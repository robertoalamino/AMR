# Preserve -c arguments before loading shared styles.
implementation=ARG1
output_file=ARG2
output_format=ARG3
grid_file=ARG4
contour_file=ARG5

# gnuplot -c plots/paper/sweep.gp IMPLEMENTATION OUTPUT FORMAT GRID CONTOURS
# GRID is ordered by alpha, with a blank line between alpha blocks.
if (!(implementation eq "corrected" || implementation eq "historical")) exit error "Invalid implementation"
image_width=1000
image_height=850
load "plots/paper/style.gp"
# Read both the tab-separated simulation grid and space-separated contours.
set datafile separator whitespace

# Contours are calculated on the recorded 11 x 11 grid.
set contour base
set cntrparam levels incremental 0,0.02,0.30
unset surface
unset key
set table contour_file
splot grid_file using 1:2:3
unset table
unset contour
set surface explicit

set title sprintf("%s sweep: paper-style view", implementation eq "corrected" ? "Corrected" : "Historical baseline")
set xlabel "{/Symbol a}"
set ylabel "{/Symbol b}"
set cblabel "Mean final q"
set xrange [-1:1]
set yrange [-1:1]
set xtics 0.5
set ytics 0.5
set cbrange [0:0.30]
set cbtics ("0" 0, "0.15" 0.15, "0.30" 0.30)
set palette defined (0 "#304fe8", 0.15 "#ffffff", 0.22 "#fff66b", 0.30 "#cf2027")
set size ratio -1
set pm3d map
set pm3d at b
set pm3d interpolate 20,20

# Keep this last: "set pm3d at b" enables implicit surfaces.
set pm3d explicit
unset grid
set label 1 "3000 steps; 100 realisations per pair; m = 0; seed 1" at screen 0.5,0.025 center font ",11"
splot grid_file using 1:2:3 with pm3d notitle, \
      contour_file using 1:2:3 with lines lc rgb "#555555" lw 0.6 notitle
unset output
