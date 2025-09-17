
set xlabel "N"
set ylabel "|error|"

set logscale x 10
set logscale y 10

plot "bad_x.dat" using 1:2 with lines title "trapez", "bad_x.dat" using 1:3 with lines title "simpson", "bad_x.dat" using 1:4 with lines title "gauss"
