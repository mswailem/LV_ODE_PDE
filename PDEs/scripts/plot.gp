set terminal pngcairo size 800,800 enhanced
unset key
system ("mkdir -p ../images")
system("rm -f ../images/*")
rgb(r,g,b) = int(r)*65536 + int(g)*256 + int(b)
datafile = '< cat *_x.dat'
stats datafile u 3 prefix 'u' nooutput
stats datafile u 4 prefix 'v' nooutput
ts=system("ls -v | grep '[[:digit:]]' | grep x | sed 's/.dat//' | sed 's/_x//'")
do for [t in ts] {
	set output '../images/'.t.'.png'
	#stats t.'_x.dat' u 3 prefix 'u' nooutput
	#stats t.'_x.dat' u 4 prefix 'v' nooutput
	set title "t=".t
	system("echo -n 'Processing time step: ".t."\r'")
	plot t.'_x.dat' u 1:2:(rgb(($3/(u_max))*255,0,($4/v_max)*255)) w lines lc rgb variable notitle
}
