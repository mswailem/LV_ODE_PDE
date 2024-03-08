@EPS
set lmargin 8
files=system("ls ../../output/chaos_diagram/ | sed 's/.dat//'")
set xlabel '{/Symbol a}'
set ylabel "fp" offset 2,0
set tmargin 3
unset key
do for [file in files] {
    set output file.'.eps'
    set title file
    plot '../../output/chaos_diagram/'.file.'.dat' notitle
}
