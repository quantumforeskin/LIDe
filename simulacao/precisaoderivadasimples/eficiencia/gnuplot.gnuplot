set terminal png size 1200,600
set output "eficiencia.png"

set xlabel "t(s)"
set ylabel "eficiencia"

set xrange [5e-6:2e-3]
set yrange [1:2]

set logscale x 10

plot "eficiencia.bonito" u 2:3
