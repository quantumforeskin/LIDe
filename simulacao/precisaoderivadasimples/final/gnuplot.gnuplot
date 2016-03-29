set terminal png size 1200,600

set xlabel "z(m)"
set ylabel "E(J)"

f(x)=2.5e5+2.29e12/sqrt(200*pi/9)*exp(-(x-25)*(x-25)/(200/9))
g(x)=-9.86664e9*x*exp(-x*x*9/200)

set output "metodo2a_final.png"
plot "evolution_final.run" u 1:3, g(x)

set ylabel "dp/dx (Pa/m)"

set output "metodo2b_final.png"
plot "evolution_final.run" u 1:($2==0?$3/g($1):1/0) w lines

