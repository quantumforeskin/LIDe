set terminal gif animate delay 10
set output "density.gif"

set xlabel "x(m)"
set ylabel "rho"

set xrange [0:2]
set yrange [200:400]

do for [i=0:999]{
  tm=i*0.01/199
  plot "evolution1.run" u 2:($1==i?$3:1/0) w lines t sprintf("%f",tm)
}

set terminal gif animate delay 10
set output "momentum.gif"

set ylabel "rhou"
set yrange [-500:500]

do for [i=0:999]{
  tm=i*0.01/199
  plot "evolution1.run" u 2:($1==i?$4:1/0) w lines t sprintf("%f",tm)
}

set terminal gif animate delay 10
set output "energy.gif"

set ylabel "Energy"
set yrange [250000:300000]

do for [i=0:999]{
  tm=i*0.01/199
  plot "evolution1.run" u 2:($1==i?$7:1/0) w lines t sprintf("%f",tm)
}

set terminal gif animate delay 10
set output "pressure.gif"

set ylabel "Pressure"
set yrange [100000:120000]

do for [i=0:999]{
  tm=i*0.01/199
  plot "evolution1.run" u 2:($1==i?$8:1/0) w lines t sprintf("%f",tm)
}




















### FOR 3D GRAPH ###

#set terminal png size 1200,600

#set xlabel "t(s)"
#set ylabel "x(m)"
#set zlabel "rho"

#set xrange [0:0.01]

#set hidden3d
#set dgrid3d 50,50 qnorm 2

#set output "rho.png"
#splot "evolution.run" u 1:2:5 w lines

#set output "rhou.png"
#splot "evolution.run" u 1:2:6 w lines

#set output "energy.png"
#splot "evolution.run" u 1:2:9 w lines

