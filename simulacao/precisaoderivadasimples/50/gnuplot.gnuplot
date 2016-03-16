set terminal png size 1200,600

set xlabel "z(m)"
set ylabel "E(J)"

f(x)=2.5e5+2.29e12/sqrt(200*pi/9)*exp(-(x-25)*(x-25)/(200/9))
g(x)=.4*2.29e12/sqrt(200*pi/9)*(25-x)/(200/9)*exp(-(x-25)*(x-25)/(200/9))

set output "initialcondition_50.png"
plot "evolution3_50.run" u 2:($1==0?$7:1/0), f(x)

set output "metodo1a_50.png"
plot "evolution3_50.run" u 2:($1==1?$4:1/0), g(x)

set output "metodo2a_50.png"
plot "evolution3_50.run" u 2:($1==2?$4:1/0), g(x)

set ylabel "dp/dx (Pa/m)"

set output "metodo1b_50.png"
plot "evolution3_50.run" u 2:($1==1?$4/g($1):1/0)

set output "metodo2b_50.png"
plot "evolution3_50.run" u 2:($1==2?$4/g($1):1/0)

