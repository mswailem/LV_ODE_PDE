@EPS
set xlabel 'n'
set ylabel '{/Symbol k}_1' offset 2,0
set output 'SH_tongue.eps'
set lmargin 8
unset key
set xtics 0.46,0.02,0.56
set xrange[0.4599:0.56]
plot '../output/fm/k1_vs_n/k0=0.245_ustar=1_vtar=1_wn=0.dat' lc rgb "red" notitle
