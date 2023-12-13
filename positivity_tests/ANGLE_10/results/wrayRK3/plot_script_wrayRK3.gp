set terminal pdfcairo enhanced color size 5,5
set output "ANGLE/plots/positivity_wrayRK3_t_37.6991.pdf"
set size square
set xlabel "t [s]"
set ylabel "y_{wrayRK3}(t)"
set xrange [0:37.6991]
set yrange [-2:2]
set grid
plot 'ANGLE/results/wrayRK3/testPositivity_wrayRK3_h_0.4_N_94.dat' using 1:2 with linespoints title 'h = 0.4',\
'ANGLE/results/wrayRK3/testPositivity_wrayRK3_h_0.8_N_47.dat' using 1:2 with linespoints title 'h = 0.8',\
'ANGLE/results/wrayRK3/testPositivity_wrayRK3_h_1.2_N_31.dat' using 1:2 with linespoints title 'h = 1.2',\
'ANGLE/results/wrayRK3/testPositivity_wrayRK3_h_2.4_N_16.dat' using 1:2 with linespoints title 'h = 2.4',\
'ANGLE/results/wrayRK3/testPositivity_wrayRK3_h_3.6_N_10.dat' using 1:2 with linespoints title 'h = 3.6',\
