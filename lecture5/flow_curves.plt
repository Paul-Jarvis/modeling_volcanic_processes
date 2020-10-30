#Gnuplot script to plot a flow curve for a Newtonian fluid

set terminal epslatex standalone

set output 'flow_curves.tex'

unset key

unset xtics
unset ytics

set xlabel '\Huge $\dot{\epsilon}$'
set ylabel '\Huge $\tau$'

set xrange[0:1]

set rmargin 16

set label '\Huge{Newtonian}' at 1, 0.97
set label '\textcolor{red}{\Huge{Bingham}}' at 1, 1.37
set label '\textcolor{red}{\Huge{$\tau_{0}$}}' at -0.09, 0.4
plot x w l lc 0 lw 3, x + 0.4 w l lc 7 dt 4 lw 3