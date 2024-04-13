@EPS
set output 'n=1.eps'
set xlabel 'a_0'
set ylabel 'b_0' offset 3
set lmargin 5
set rmargin 0.8
set yrange[-0.001:0.401]
set ytics 0,0.1,0.4
unset key
plot '../../output/fm/ustar_vs_vstar/n=1.dat' u 1:2 lc rgb "red" notitle
