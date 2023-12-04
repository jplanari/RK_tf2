set terminal pdf
set output "plots/positivity_heunRK2_t_10.pdf"
set size square
set xlabel "t [s]"
set ylabel "y_{heunRK2}(t)"
set xrange [0:10]
set grid
set yrange [-3:4]
plot 'results/heunRK2/testPositivity_heunRK2_h_0.1_N_100.dat' using 1:2 title 'h = 0.1',\
'results/heunRK2/testPositivity_heunRK2_h_0.2_N_50.dat' using 1:2 title 'h = 0.2',\
'results/heunRK2/testPositivity_heunRK2_h_0.4_N_25.dat' using 1:2 title 'h = 0.4',\
'results/heunRK2/testPositivity_heunRK2_h_0.8_N_12.dat' using 1:2 title 'h = 0.8',\
'results/heunRK2/testPositivity_heunRK2_h_1.6_N_6.dat' using 1:2 title 'h = 1.6',\
'results/heunRK2/testPositivity_heunRK2_h_2.4_N_4.dat' using 1:2 title 'h = 2.4',\
'results/heunRK2/testPositivity_heunRK2_h_3.2_N_3.dat' using 1:2 title 'h = 3.2',\
'results/heunRK2/testPositivity_heunRK2_h_4.0_N_2.dat' using 1:2 title 'h = 4.0',\
