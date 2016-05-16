
#4 d density
#5 s speed
#8 e energy
#9 p pressure

set xlabel "x(m)"

n=10

set terminal gif animate delay 100
set ylabel "p(Kgm^{-3})"
set yrange [0:5]
set output "10d.gif"
i=0
do for[i=0:n-1:1]{
plot "10.run" u 3:($1==i?$4:1/0) w lines title sprintf("t=%.5es",1.11111e-5*i)
}

set terminal gif animate delay 100
set ylabel "v(ms^{-1})"
set yrange [-1e6:1e6]
set output "10v.gif"
i=0
do for[i=0:n-1:1]{
plot "10.run" u 3:($1==i?$5:1/0) w lines title sprintf("t=%.5es",1.11111e-5*i)
}

set terminal gif animate delay 100
set ylabel "E(J)"
set yrange [0:5e11]
set output "10e.gif"
i=0
do for[i=0:n-1:1]{
plot "10.run" u 3:($1==i?$8:1/0) w lines title sprintf("t=%.5es",1.11111e-5*i)
}

set terminal gif animate delay 100
set ylabel "P(Pa)"
set yrange [1e3:1.5e11]
set output "10p.gif"
i=0
do for[i=0:n-1:1]{
plot "10.run" u 3:($1==i?$9:1/0) w lines title sprintf("t=%.5es",1.11111e-5*i)
}

n=20

set terminal gif animate delay 100
set ylabel "p(Kgm^{-3})"
set yrange [0:5]
set output "20d.gif"
i=0
do for[i=0:n-1:1]{
plot "20.run" u 3:($1==i?$4:1/0) w lines title sprintf("t=%.5es",5.26316e-6*i)
}

set terminal gif animate delay 100
set ylabel "v(ms^{-1})"
set yrange [-1e6:1e6]
set output "20v.gif"
i=0
do for[i=0:n-1:1]{
plot "20.run" u 3:($1==i?$5:1/0) w lines title sprintf("t=%.5es",5.26316e-6*i)
}

set terminal gif animate delay 100
set ylabel "E(J)"
set yrange [0:5e11]
set output "20e.gif"
i=0
do for[i=0:n-1:1]{
plot "20.run" u 3:($1==i?$8:1/0) w lines title sprintf("t=%.5es",5.26316e-6*i)
}

set terminal gif animate delay 100
set ylabel "P(Pa)"
set yrange [1e3:1.5e11]
set output "20p.gif"
i=0
do for[i=0:n-1:1]{
plot "20.run" u 3:($1==i?$9:1/0) w lines title sprintf("t=%.5es",5.26316e-6*i)
}

n=30

set terminal gif animate delay 100
set ylabel "p(Kgm^{-3})"
set yrange [0:5]
set output "30d.gif"
i=0
do for[i=0:n-1:1]{
plot "30.run" u 3:($1==i?$4:1/0) w lines title sprintf("t=%.5es",3.44828e-6*i)
}

set terminal gif animate delay 100
set ylabel "v(ms^{-1})"
set yrange [-1e6:1e6]
set output "30v.gif"
i=0
do for[i=0:n-1:1]{
plot "30.run" u 3:($1==i?$5:1/0) w lines title sprintf("t=%.5es",3.44828e-6*i)
}

set terminal gif animate delay 100
set ylabel "E(J)"
set yrange [0:5e11]
set output "30e.gif"
i=0
do for[i=0:n-1:1]{
plot "30.run" u 3:($1==i?$8:1/0) w lines title sprintf("t=%.5es",3.44828e-6*i)
}

set terminal gif animate delay 100
set ylabel "P(Pa)"
set yrange [1e3:1.5e11]
set output "30p.gif"
i=0
do for[i=0:n-1:1]{
plot "30.run" u 3:($1==i?$9:1/0) w lines title sprintf("t=%.5es",3.44828e-6*i)
}

n=50

set terminal gif animate delay 100
set ylabel "p(Kgm^{-3})"
set yrange [0:5]
set output "50d.gif"
i=0
do for[i=0:n-1:1]{
plot "50.run" u 3:($1==i?$4:1/0) w lines title sprintf("t=%.5es",2.04082e-06*i)
}

set terminal gif animate delay 100
set ylabel "v(ms^{-1})"
set yrange [-1e6:1e6]
set output "50v.gif"
i=0
do for[i=0:n-1:1]{
plot "50.run" u 3:($1==i?$5:1/0) w lines title sprintf("t=%.5es",2.04082e-06*i)
}

set terminal gif animate delay 100
set ylabel "E(J)"
set yrange [0:5e11]
set output "50e.gif"
i=0
do for[i=0:n-1:1]{
plot "50.run" u 3:($1==i?$8:1/0) w lines title sprintf("t=%.5es",2.04082e-06*i)
}

set terminal gif animate delay 100
set ylabel "P(Pa)"
set yrange [1e3:1.5e11]
set output "50p.gif"
i=0
do for[i=0:n-1:1]{
plot "50.run" u 3:($1==i?$9:1/0) w lines title sprintf("t=%.5es",2.04082e-06*i)
}

n=75

set terminal gif animate delay 100
set ylabel "p(Kgm^{-3})"
set yrange [0:5]
set output "75d.gif"
i=0
do for[i=0:n-1:1]{
plot "75.run" u 3:($1==i?$4:1/0) w lines title sprintf("t=%.5es",1.35135e-06*i)
}

set terminal gif animate delay 100
set ylabel "v(ms^{-1})"
set yrange [-1e6:1e6]
set output "75v.gif"
i=0
do for[i=0:n-1:1]{
plot "75.run" u 3:($1==i?$5:1/0) w lines title sprintf("t=%.5es",1.35135e-06*i)
}

set terminal gif animate delay 100
set ylabel "E(J)"
set yrange [0:5e11]
set output "75e.gif"
i=0
do for[i=0:n-1:1]{
plot "75.run" u 3:($1==i?$8:1/0) w lines title sprintf("t=%.5es",1.35135e-06*i)
}

set terminal gif animate delay 100
set ylabel "P(Pa)"
set yrange [1e3:1.5e11]
set output "75p.gif"
i=0
do for[i=0:n-1:1]{
plot "75.run" u 3:($1==i?$9:1/0) w lines title sprintf("t=%.5es",1.35135e-06*i)
}

n=50#100

set terminal gif animate delay 100
set ylabel "p(Kgm^{-3})"
set yrange [0:5]
set output "100d.gif"
i=0
do for[i=0:n-1:1]{
plot "100.run" u 3:($1==2*i?$4:1/0) w lines title sprintf("t=%.5es",2.0202e-06*i)
}

set terminal gif animate delay 100
set ylabel "v(ms^{-1})"
set yrange [-1e6:1e6]
set output "100v.gif"
i=0
do for[i=0:n-1:1]{
plot "100.run" u 3:($1==2*i?$5:1/0) w lines title sprintf("t=%.5es",2.0202e-06*i)
}

set terminal gif animate delay 100
set ylabel "E(J)"
set yrange [0:5e11]
set output "100e.gif"
i=0
do for[i=0:n-1:1]{
plot "100.run" u 3:($1==2*i?$8:1/0) w lines title sprintf("t=%.5es",2.0202e-06*i)
}

set terminal gif animate delay 100
set ylabel "P(Pa)"
set yrange [1e3:1.5e11]
set output "100p.gif"
i=0
do for[i=0:n-1:1]{
plot "100.run" u 3:($1==2*i?$9:1/0) w lines title sprintf("t=%.5es",2.0202e-06*i)
}

n=50#200

set terminal gif animate delay 100
set ylabel "p(Kgm^{-3})"
set yrange [0:5]
set output "200d.gif"
i=0
do for[i=0:n-1:1]{
plot "200.run" u 3:($1==4*i?$4:1/0) w lines title sprintf("t=%.5es",2.01005e-06*i)
}

set terminal gif animate delay 100
set ylabel "v(ms^{-1})"
set yrange [-1e6:1e6]
set output "200v.gif"
i=0
do for[i=0:n-1:1]{
plot "200.run" u 3:($1==4*i?$5:1/0) w lines title sprintf("t=%.5es",2.01005e-06*i)
}

set terminal gif animate delay 100
set ylabel "E(J)"
set yrange [0:5e11]
set output "200e.gif"
i=0
do for[i=0:n-1:1]{
plot "200.run" u 3:($1==4*i?$8:1/0) w lines title sprintf("t=%.5es",2.01005e-06*i)
}

set terminal gif animate delay 100
set ylabel "P(Pa)"
set yrange [1e3:1.5e11]
set output "200p.gif"
i=0
do for[i=0:n-1:1]{
plot "200.run" u 3:($1==4*i?$9:1/0) w lines title sprintf("t=%.5es",2.01005e-06*i)
}

n=50#300

set terminal gif animate delay 100
set ylabel "p(Kgm^{-3})"
set yrange [0:5]
set output "300d.gif"
i=0
do for[i=0:n-1:1]{
plot "300.run" u 3:($1==6*i?$4:1/0) w lines title sprintf("t=%.5es",2.00669e-06*i)
}

set terminal gif animate delay 100
set ylabel "v(ms^{-1})"
set yrange [-1e6:1e6]
set output "300v.gif"
i=0
do for[i=0:n-1:1]{
plot "300.run" u 3:($1==6*i?$5:1/0) w lines title sprintf("t=%.5es",2.00669e-06*i)
}

set terminal gif animate delay 100
set ylabel "E(J)"
set yrange [0:5e11]
set output "300e.gif"
i=0
do for[i=0:n-1:1]{
plot "300.run" u 3:($1==6*i?$8:1/0) w lines title sprintf("t=%.5es",2.00669e-06*i)
}

set terminal gif animate delay 100
set ylabel "P(Pa)"
set yrange [1e3:1.5e11]
set output "300p.gif"
i=0
do for[i=0:n-1:1]{
plot "300.run" u 3:($1==6*i?$9:1/0) w lines title sprintf("t=%.5es",2.00669e-06*i)
}

n=50#0

set terminal gif animate delay 100
set ylabel "p(Kgm^{-3})"
set yrange [0:5]
set output "500d.gif"
i=0
do for[i=0:n-1:1]{
plot "500.run" u 3:($1==10*i?$4:1/0) w lines title sprintf("t=%.5es",2.00401e-06*i)
}

set terminal gif animate delay 100
set ylabel "v(ms^{-1})"
set yrange [-1e6:1e6]
set output "500v.gif"
i=0
do for[i=0:n-1:1]{
plot "500.run" u 3:($1==10*i?$5:1/0) w lines title sprintf("t=%.5es",2.00401e-06*i)
}

set terminal gif animate delay 100
set ylabel "E(J)"
set yrange [0:5e11]
set output "500e.gif"
i=0
do for[i=0:n-1:1]{
plot "500.run" u 3:($1==10*i?$8:1/0) w lines title sprintf("t=%.5es",2.00401e-06*i)
}

set terminal gif animate delay 100
set ylabel "P(Pa)"
set yrange [1e3:1.5e11]
set output "500p.gif"
i=0
do for[i=0:n-1:1]{
plot "500.run" u 3:($1==10*i?$9:1/0) w lines title sprintf("t=%.5es",2.00401e-06*i)
}

n=50#750

set terminal gif animate delay 100
set ylabel "p(Kgm^{-3})"
set yrange [0:5]
set output "750d.gif"
i=0
do for[i=0:n-1:1]{
plot "750.run" u 3:($1==15*i?$4:1/0) w lines title sprintf("t=%.5es",2.00267e-06*i)
}

set terminal gif animate delay 100
set ylabel "v(ms^{-1})"
set yrange [-1e6:1e6]
set output "750v.gif"
i=0
do for[i=0:n-1:1]{
plot "750.run" u 3:($1==15*i?$5:1/0) w lines title sprintf("t=%.5es",2.00267e-06*i)
}

set terminal gif animate delay 100
set ylabel "E(J)"
set yrange [0:5e11]
set output "750e.gif"
i=0
do for[i=0:n-1:1]{
plot "750.run" u 3:($1==15*i?$8:1/0) w lines title sprintf("t=%.5es",2.00267e-06*i)
}

set terminal gif animate delay 100
set ylabel "P(Pa)"
set yrange [1e3:1.5e11]
set output "750p.gif"
i=0
do for[i=0:n-1:1]{
plot "750.run" u 3:($1==15*i?$9:1/0) w lines title sprintf("t=%.5es",2.00267e-06*i)
}

n=50#1000

set terminal gif animate delay 100
set ylabel "p(Kgm^{-3})"
set yrange [0:5]
set output "1000d.gif"
i=0
do for[i=0:n-1:1]{
plot "1000.run" u 3:($1==20*i?$4:1/0) w lines title sprintf("t=%.5es",2.002e-06*i)
}

set terminal gif animate delay 100
set ylabel "v(ms^{-1})"
set yrange [-1e6:1e6]
set output "1000v.gif"
i=0
do for[i=0:n-1:1]{
plot "1000.run" u 3:($1==20*i?$5:1/0) w lines title sprintf("t=%.5es",2.002e-06*i)
}

set terminal gif animate delay 100
set ylabel "E(J)"
set yrange [0:5e11]
set output "1000e.gif"
i=0
do for[i=0:n-1:1]{
plot "1000.run" u 3:($1==20*i?$8:1/0) w lines title sprintf("t=%.5es",2.002e-06*i)
}

set terminal gif animate delay 100
set ylabel "P(Pa)"
set yrange [1e3:1.5e11]
set output "1000p.gif"
i=0
do for[i=0:n-1:1]{
plot "1000.run" u 3:($1==20*i?$9:1/0) w lines title sprintf("t=%.5es",2.002e-06*i)
}

n=50#2000

set terminal gif animate delay 100
set ylabel "p(Kgm^{-3})"
set yrange [0:5]
set output "2000d.gif"
i=0
do for[i=0:n-1:1]{
plot "2000.run" u 3:($1==40*i?$4:1/0) w lines title sprintf("t=%.5es",2.001e-06*i)
}

set terminal gif animate delay 100
set ylabel "v(ms^{-1})"
set yrange [-1e6:1e6]
set output "2000v.gif"
i=0
do for[i=0:n-1:1]{
plot "2000.run" u 3:($1==40*i?$5:1/0) w lines title sprintf("t=%.5es",2.001e-06*i)
}

set terminal gif animate delay 100
set ylabel "E(J)"
set yrange [0:5e11]
set output "2000e.gif"
i=0
do for[i=0:n-1:1]{
plot "2000.run" u 3:($1==40*i?$8:1/0) w lines title sprintf("t=%.5es",2.001e-06*i)
}

set terminal gif animate delay 100
set ylabel "P(Pa)"
set yrange [1e3:1.5e11]
set output "2000p.gif"
i=0
do for[i=0:n-1:1]{
plot "2000.run" u 3:($1==40*i?$9:1/0) w lines title sprintf("t=%.5es",2.001e-06*i)
}

n=50#3000

set terminal gif animate delay 100
set ylabel "p(Kgm^{-3})"
set yrange [0:5]
set output "3000d.gif"
i=0
do for[i=0:n-1:1]{
plot "3000.run" u 3:($1==60*i?$4:1/0) w lines title sprintf("t=%.5es",2.00067e-06*i)
}

set terminal gif animate delay 100
set ylabel "v(ms^{-1})"
set yrange [-1e6:1e6]
set output "3000v.gif"
i=0
do for[i=0:n-1:1]{
plot "3000.run" u 3:($1==60*i?$5:1/0) w lines title sprintf("t=%.5es",2.00067e-06*i)
}

set terminal gif animate delay 100
set ylabel "E(J)"
set yrange [0:5e11]
set output "3000e.gif"
i=0
do for[i=0:n-1:1]{
plot "3000.run" u 3:($1==60*i?$8:1/0) w lines title sprintf("t=%.5es",2.00067e-06*i)
}

set terminal gif animate delay 100
set ylabel "P(Pa)"
set yrange [1e3:1.5e11]
set output "3000p.gif"
i=0
do for[i=0:n-1:1]{
plot "3000.run" u 3:($1==60*i?$9:1/0) w lines title sprintf("t=%.5es",2.00067e-06*i)
}


