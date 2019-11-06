#Gnuplot script to plot a flow curve for a Newtonian fluid

set terminal epslatex standalone

set output 'flow_curves.tex'

unset key

unset xtics
unset ytics

set xlabel '\huge $\dot{\epsilon}$'
set ylabel '\huge $\tau$'

set title '\huge $\tau = \eta \dot{\epsilon}$'

set xrange[0:1]

set rmargin 14

set label '\textcolor{blue}{\large{Shear thickening}}' at 1, 0.85
set label '\large{Newtonian}' at 1, 0.97
set label '\textcolor{green}{\large{Shear thinning}}' at 1, 1.07
set label '\textcolor{red}{\large{Bingham}}' at 1, 1.37
set label '\textcolor{red}{\large{$\tau_{0}$}}' at -0.05, 0.4
plot x w l lc 0 lw 3, 0.9 * x**2 w l lc 6 dt 2 lw 3, 1.1 * x**(1.0/2.0) w l lc 2 dt 3 lw 3, x + 0.4 w l lc 7 dt 4 lw 3