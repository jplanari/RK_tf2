set terminal pdfcairo enhanced color size 5,5
set output "ANGLE/plots/positivity_paramEuler_t_37.6991.pdf"
set size square
set xlabel "t [s]"
set ylabel "y_{paramEuler}(t)"
set xrange [0:37.6991]
set yrange [-2:2]
set grid
plot 'ANGLE/results/paramEuler/testPositivity_paramEuler_h_0.8_N_47.dat' using 1:2 with linespoints title 'h = 0.8',\
'ANGLE/results/paramEuler/testPositivity_paramEuler_h_1.8_N_21.dat' using 1:2 with linespoints title 'h = 1.8',\
'ANGLE/results/paramEuler/testPositivity_paramEuler_h_2.3_N_16.dat' using 1:2 with linespoints title 'h = 2.3',\
'ANGLE/results/paramEuler/testPositivity_paramEuler_h_2.4_N_16.dat' using 1:2 with linespoints title 'h = 2.4',\
'ANGLE/results/paramEuler/testPositivity_paramEuler_h_2.5_N_15.dat' using 1:2 with linespoints title 'h = 2.5',\
'ANGLE/results/paramEuler/testPositivity_paramEuler_h_2.75_N_14.dat' using 1:2 with linespoints title 'h = 2.75',\
