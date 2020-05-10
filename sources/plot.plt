# gnuplot -c plot.plt -p
set title "Lotka-Volterra Equations"
set grid
set xlabel "Time"
set ylabel "Population"
plot "output.txt" using 1:2 with lines title "Prey", "output.txt" using 1:3 with lines title "Predator"
