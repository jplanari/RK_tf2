set terminal pdfcairo enhanced color size 5,5
set output "DIFFUSIVE/plots/positivity_stdRK4_t_15.pdf"
set size square
set xlabel "t [s]" font ',25'
set ylabel "y_{stdRK4}(t)" font ',25'
set xrange [0:15]
set yrange [-1:1]
set key font ',25'

set xtics font ',20'
set ytics font ',20'

set grid
plot 'DIFFUSIVE/results/stdRK4/testPositivity_stdRK4_h_0.8_N_19.dat' using 1:2 with linespoints title 'h = 0.8',\
'DIFFUSIVE/results/stdRK4/testPositivity_stdRK4_h_1.2_N_12.dat' using 1:2 with linespoints title 'h = 1.2',\
'DIFFUSIVE/results/stdRK4/testPositivity_stdRK4_h_1.8_N_8.dat' using 1:2 with linespoints title 'h = 1.8',\
'DIFFUSIVE/results/stdRK4/testPositivity_stdRK4_h_2.3_N_6.dat' using 1:2 with linespoints title 'h = 2.3',\
'DIFFUSIVE/results/stdRK4/testPositivity_stdRK4_h_2.6_N_6.dat' using 1:2 with linespoints title 'h = 2.6',\
