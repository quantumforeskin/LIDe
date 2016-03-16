set terminal png size 1200,600
set output "initial condition.png"

set xlabel "z(m)"
set ylabel "E(J)"

o=10./3
f(x)=2.29e12/sqrt(2*pi*o*o)*exp(-(x-29700)*(x-29700)/(2*o*o))

plot "evolution3.run" u 2:7 w lines, f(x)
