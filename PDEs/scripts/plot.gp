set terminal png
unset key
system ("mkdir -p ../images")
system("rm ../images/*")
rgb(r,g,b) = int(r)*65536 + int(g)*256 + int(b)
ts=system("ls -v | grep '[[:digit:]]' | grep x | sed 's/.dat//' | sed 's/_x//'")
do for [t in ts] {
	set output '../images/'.t.'.png'
	set title "t=".t
	plot t.'_x.dat' u 1:2:(rgb($3*255,0,$4*255)) pt 7 lc rgb variable notitle
}
