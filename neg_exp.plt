
set xlabel "N"
set ylabel "|error|"

set logscale x 10
set logscale y 10

plot "neg_exp.dat" using 1:2 with lines title "trapez", "neg_exp.dat" using 1:3 with lines title "simpson", "neg_exp.dat" using 1:4 with lines title "gauss"
