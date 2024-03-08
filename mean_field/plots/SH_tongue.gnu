@EPS
set xlabel 'n'
set ylabel 'k_1' offset 2,0
set output 'SH_tongue.eps'
set lmargin 8
unset key
plot '../output/fm/a0=1_b0=1_k0=0.246.dat' notitle
