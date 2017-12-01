# Gnuplot script file for plotting data in file
set terminal png size 400,300 enhanced font "Helvetica,20"
set output 'output.png'

set xtics font ",20"
set ytics font ",20"

set xrange [-1.0:20]
set yrange [12:16]

set key font ",20" spacing 1.5
set lmargin at screen 0.15
set rmargin at screen 0.95
plot "300K" u ($2/1000):($11) w l lw 2 title "300K", \
     "1000K" u ($2/1000):($11) w l lw 2 title "1000K", \
     "2000K" u ($2/1000):($11) w l lw 2 title "2000K", \
     "4000K" u ($2/1000):($11) w l lw 2 title "4000K", \
     "6000K" u ($2/1000):($11) w l lw 2 title "6000K"
