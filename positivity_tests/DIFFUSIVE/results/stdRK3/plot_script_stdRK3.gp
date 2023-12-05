set terminal pdfcairo enhanced color size 5,5
set output "DIFFUSIVE/plots/positivity_stdRK3_t_37.6991.pdf"
set size square
set xlabel "t [s]"
set ylabel "y_{stdRK3}(t)"
set xrange [0:37.6991]
set yrange [-2:2]
set grid
plot 'DIFFUSIVE/results/stdRK3/testPositivity_stdRK3_h_0.8_N_47.dat' using 1:2 with linespoints title 'h = 0.8',\
'DIFFUSIVE/results/stdRK3/testPositivity_stdRK3_h_1.8_N_21.dat' using 1:2 with linespoints title 'h = 1.8',\
'DIFFUSIVE/results/stdRK3/testPositivity_stdRK3_h_2.3_N_16.dat' using 1:2 with linespoints title 'h = 2.3',\
'DIFFUSIVE/results/stdRK3/testPositivity_stdRK3_h_2.4_N_16.dat' using 1:2 with linespoints title 'h = 2.4',\
'DIFFUSIVE/results/stdRK3/testPositivity_stdRK3_h_2.5_N_15.dat' using 1:2 with linespoints title 'h = 2.5',\
'DIFFUSIVE/results/stdRK3/testPositivity_stdRK3_h_2.75_N_14.dat' using 1:2 with linespoints title 'h = 2.75',\
