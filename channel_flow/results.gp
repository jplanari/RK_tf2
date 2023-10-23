
set terminal eps enhanced color
set output "cf_180_64.eps"

# Set your plot title and labels

set xlabel "y^+"
set ylabel "u^+(y^+)"

set xrange [0.5:180]
set yrange [0.5:50]

set grid
set logscale xy

set multiplot

plot "results/profile_noModel_cells_64.dat" using 2:3 title "Re = 180"

replot "results/profile_noModel_cells_64.dat" using 2:2 with lines title "u^+=y^"

replot "results/profile_noModel_cells_64.dat" using 2:(5.5+2.5*log($2)) with lines title "Log law"

unset multiplot
