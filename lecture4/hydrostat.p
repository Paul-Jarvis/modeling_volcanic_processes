#Gnuplot script to plot two hydrostatic pressure profiles

set terminal epslatex standalone

set output 'hydrostat.tex'

set yrange[1:0]
set xrange[0:1]

set border 0
unset xtics
unset ytics
set arrow from 0,0 to 1, 0 lw 5
set arrow from 0,0 to 0, 1 lw 5 

set key spacing 2

set x2label '\Huge $P$'
set ylabel '\Huge $z$'

plot x title '\Huge $P_{2}$' lc 6  lw 5, 3*x/4 title '\Huge $P_{1}$' lc 7  lw 5

set output 'hydrostat_diff.tex'

plot x title '\Huge $P_{2}$' lc 6  lw 5, 3*x/4 title '\Huge $P_{1}$' lc 7  lw 5, 4 * x title '\Huge $\Delta P$' lc 0 lw 5

set output 'hydrostat_1fluid.tex'

plot x notitle lc 6 lw 5
