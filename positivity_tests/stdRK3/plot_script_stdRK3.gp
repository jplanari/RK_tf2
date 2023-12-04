set terminal pdf
set output "positivity_stdRK3_t_5.0.pdf"
set xlabel "t [s]"
set ylabel "y_{stdRK3}(t)"
set xrange [0:5.0]
set grid

plot \

'testPositivity_stdRK3_h_0.1_N_50.dat' using 1:2 title 'h = 0.1'
,\
'testPositivity_stdRK3_h_0.2_N_25.dat' using 1:2 title 'h = 0.2'
,\
'testPositivity_stdRK3_h_0.4_N_12.dat' using 1:2 title 'h = 0.4'
,\
'testPositivity_stdRK3_h_0.8_N_6.dat' using 1:2 title 'h = 0.8'
,\
'testPositivity_stdRK3_h_1.6_N_3.dat' using 1:2 title 'h = 1.6'
,\
'testPositivity_stdRK3_h_2.4_N_2.dat' using 1:2 title 'h = 2.4'
,\
'testPositivity_stdRK3_h_3.2_N_2.dat' using 1:2 title 'h = 3.2'
,\
'testPositivity_stdRK3_h_4.0_N_1.dat' using 1:2 title 'h = 4.0'
