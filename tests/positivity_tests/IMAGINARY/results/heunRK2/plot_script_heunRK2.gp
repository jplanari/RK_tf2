set terminal pdfcairo enhanced color size 5,5
set output "IMAGINARY/plots/positivity_heunRK2_t_37.6991.pdf"
set size square
set xlabel "t [s]"
set ylabel "y_{heunRK2}(t)"
set xrange [0:37.6991]
set yrange [-2:2]
set grid
plot 'IMAGINARY/results/heunRK2/testPositivity_heunRK2_h_0.8_N_47.dat' using 1:2 with linespoints title 'h = 0.8',\
'IMAGINARY/results/heunRK2/testPositivity_heunRK2_h_1.8_N_21.dat' using 1:2 with linespoints title 'h = 1.8',\
'IMAGINARY/results/heunRK2/testPositivity_heunRK2_h_2.3_N_16.dat' using 1:2 with linespoints title 'h = 2.3',\
'IMAGINARY/results/heunRK2/testPositivity_heunRK2_h_2.4_N_16.dat' using 1:2 with linespoints title 'h = 2.4',\
'IMAGINARY/results/heunRK2/testPositivity_heunRK2_h_2.5_N_15.dat' using 1:2 with linespoints title 'h = 2.5',\
'IMAGINARY/results/heunRK2/testPositivity_heunRK2_h_2.75_N_14.dat' using 1:2 with linespoints title 'h = 2.75',\
