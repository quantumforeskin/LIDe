set terminal png size 1200,600

set xlabel "z(m)"
set ylabel "E(J)"

f(x)=2.5e5+2.29e12/sqrt(200*pi/9)*exp(-(x-25)*(x-25)/(200/9))
g(x)=-9.86664e9*(x-25)*exp(-(x-25)*(x-25)*9/200)

set output "initialcondition_1000.png"
plot "evolution3_1000.run" u 2:($1==0?$7:1/0), f(x)

set output "metodo1a_1000.png"
plot "evolution3_1000.run" u 2:($1==1?$4:1/0), g(x)

set output "metodo2a_1000.png"
plot "evolution3_1000.run" u 2:($1==2?$4:1/0), g(x)

set ylabel "dp/dx (Pa/m)"
set yrange [0.5:1.5]

set output "metodo1b_1000.png"
plot "evolution3_1000.run" u 2:($1==1?$4/g($2):1/0)

set yrange [0.999:1.005]

set output "metodo2b_1000.png"
plot "evolution3_1000.run" u 2:($1==2?$4/g($2):1/0) w lines

