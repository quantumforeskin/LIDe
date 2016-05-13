set terminal gif animate
set output "caralho.gif"
n=10

set xrange [29680:29720]
set yrange [1.26:1.28]

i=0
do for[i=0:n:1]{
plot "evolution.run" u 3:($1==i?$4:1/0) w lines
}

