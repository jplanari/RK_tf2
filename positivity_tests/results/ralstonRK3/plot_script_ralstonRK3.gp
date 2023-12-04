set terminal pdf
set output "plots/positivity_ralstonRK3_t_10.pdf"
set size square
set xlabel "t [s]"
set ylabel "y_{ralstonRK3}(t)"
set xrange [0:10]
set grid
set yrange [-3:4]
plot 'results/ralstonRK3/testPositivity_ralstonRK3_h_0.1_N_100.dat' using 1:2 title 'h = 0.1',\
'results/ralstonRK3/testPositivity_ralstonRK3_h_0.2_N_50.dat' using 1:2 title 'h = 0.2',\
'results/ralstonRK3/testPositivity_ralstonRK3_h_0.4_N_25.dat' using 1:2 title 'h = 0.4',\
'results/ralstonRK3/testPositivity_ralstonRK3_h_0.8_N_12.dat' using 1:2 title 'h = 0.8',\
'results/ralstonRK3/testPositivity_ralstonRK3_h_1.6_N_6.dat' using 1:2 title 'h = 1.6',\
'results/ralstonRK3/testPositivity_ralstonRK3_h_2.4_N_4.dat' using 1:2 title 'h = 2.4',\
'results/ralstonRK3/testPositivity_ralstonRK3_h_3.2_N_3.dat' using 1:2 title 'h = 3.2',\
'results/ralstonRK3/testPositivity_ralstonRK3_h_4.0_N_2.dat' using 1:2 title 'h = 4.0',\
