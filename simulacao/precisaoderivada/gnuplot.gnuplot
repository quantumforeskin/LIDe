set terminal png size 1200,600

set xlabel "z(m)"
set ylabel "dE/dx(J/m)"

o=10./3
f(x)=-.4*2.29e12/sqrt(2*pi*10/3*10/3)*(x-29700)/(10/3*10/3)*exp(-(x-29700)*(x-29700)/(2*10/3*10/3))

set output "metodo1a.png"
plot "evolution3.run" u 2:($1==1?4:1/0) w lines, f(x)

set output "metodo2a.png"
plot "evolution3.run" u 2:($1==0?4:1/0) w lines, f(x)

set output "metodo1b.png"
plot "evolution3.run" u 2:($1==1?4/f($2):1/0) w lines

set output "metodo2b.png"
plot "evolution3.run" u 2:($1==1?4/f($2):1/0) w lines


