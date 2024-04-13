@EPS
set output 'below_tongue.eps'
set xlabel '{/Symbol a}'
set ylabel 'a(mT)' offset 2 rotate by 90
unset key
set rmargin 0.8
set lmargin 6
plot '../../output/bifurcation/alpha/for_thesis/below_SH_tongue_1.dat' u 1:2 lc rgb "blue" notitle
