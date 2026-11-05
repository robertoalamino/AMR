# Shared paper-style rendering. The calling script sets width and height.
# The caller captures output_file and output_format before loading this file.
if (output_format eq "png") {
    set terminal pngcairo size image_width,image_height enhanced font "serif,14"
} else {
    if (output_format eq "svg") {
        set terminal svg size image_width,image_height enhanced font "serif,14"
    } else {
        print "Expected png or svg output format"
        exit error "Invalid output format"
    }
}
set output output_file
set encoding utf8
set datafile separator "\t"
set datafile missing "NaN"
set border linewidth 1 lc rgb "#555555"
set tics out nomirror
set grid back lc rgb "#e5e5e5" linewidth 0.5
set xlabel "q"
set ylabel "Fraction of all realisations"
