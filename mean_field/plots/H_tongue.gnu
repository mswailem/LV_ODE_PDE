@EPS
set xlabel 'n'
set ylabel '{/Symbol k}_1' offset 2,0
set output 'H_tongue.eps'
set lmargin 8
unset key
plot '../output/fm/k1_vs_n/k0=0.095_ustar=5_vtar=0.18_wn=0.dat' lc rgb "red" notitle
