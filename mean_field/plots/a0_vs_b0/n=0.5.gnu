@EPS
set output 'n=0.5.eps'
set xlabel 'a_0'
set ylabel 'b_0' offset 2
set lmargin 5
set rmargin 0.8
set yrange[-0.001:5.001]
set ytics 0,1,5
plot '../../output/fm/ustar_vs_vstar/n=0.5.dat' u 1:2 lc rgb "red" notitle
