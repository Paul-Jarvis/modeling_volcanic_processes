#Gnuplot script to plot a flow curve for a Newtonian fluid

set terminal epslatex standalone

set output 'newtonian.tex'

unset key

unset xtics
unset ytics

set xlabel '\huge $\dot{\epsilon}$'
set ylabel '\huge $\tau$'

set title '\huge $\tau = \eta \dot{\epsilon}$'

set label '\huge Gradient = $\eta$' at 2, 1
plot x w l lc 0 