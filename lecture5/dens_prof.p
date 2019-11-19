#Gnuplot script to plot an example of the density profile of a spreading cloud

set terminal epslatex standalone

set output 'dens_prof.tex'

set xlabel 'Density'

set ylabel 'Altitude'

set xrange[2:8]

set samples 5000

f1(x) = (x < 4) ? -x:1/0
f2(x) = (x > 6) ? -x:1/0
f3(x) = (x > 4 && x < 6) ? -x:1/0
f4(x) = (x > 4 && x < 5) ? -4:1/0
f5(x) = (x > 5 && x < 6) ? -6:1/0

set arrow from 5, -6 to 5, -4 nohead lc 0 lw 2
#set arrow from 5, -6 to 5, -4 nohead lc 7 lw 3 dt 3

unset key

unset xtics

unset ytics

unset border

set arrow from 2, -8 to 8.2, -8
set arrow from 2, -8 to 2, -1.75

#set label 'Experimental density profile' at 2.5, -7
#set arrow from 5.05, -7 to 5.95, -7

set label 'density profile' at 6.5, -6
set label 'Atmospheric' at 6.55, -5.7
set arrow from 7.1, -6.1 to 7.1, -7.0
   
set label 'Mixed volcanic cloud' at 6, -4.5
set arrow from 5.95, -4.5 to 5.1, -4.5

set arrow from 5, -6.1 to 5, -6.4 nohead
set arrow from 6, -6.1 to 6, -6.4 nohead
set arrow from 5, -6.25 to 6, -6.25 heads

set label '$ \Delta \rho $' at 5.4, -6.4

#set arrow from 6, -6 to 6, -8 dt 3 lc 7 lw 3 nohead
   
plot f1(x) w l lc 0 lw 2,\
     f2(x) w l lc 0 lw 2,\
     f3(x) w l lc 0 dt 2 lw 2,\
     f4(x) w l lc 0 lw 2,\
     f5(x) w l lc 0 lw 2
