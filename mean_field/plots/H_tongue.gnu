@EPS
set xlabel 'n'
set ylabel 'k_1' offset 2,0
set output 'H_tongue.eps'
set lmargin 8
unset key
plot '../output/fm/a0=5_b0=0.18_k0=0.095.dat' notitle
