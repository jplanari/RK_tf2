
set terminal eps enhanced color
set output "cf_180_64.eps"

# Set your plot title and labels

set xlabel "y^+"
set ylabel "u^+(y^+)"

set xrange [0.5:180]
set yrange [0.5:50]

set grid
set logscale x

set multiplot

plot ARG1 using 2:3 title "Re = 180"
replot x with lines title "u^+=y^+"
replot 5.5+2.5*log(x) with lines title "Log law"

unset multiplot
