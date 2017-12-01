# Gnuplot script file for plotting data in file
set terminal png size 400,300 enhanced font "Helvetica,20"
set output 'Au147.png'

set xtics font ",20"
set ytics font ",20"

set xrange [1:6]
set yrange [0:0.6]

set key font ",20" spacing 1.5
set lmargin at screen 0.15
set rmargin at screen 0.95
plot "Au147_DFT.dat" u 1:2 w lp pt 5 ps 1 title "Au147_DFT", \
     "sc_exafs.chir" u 1:2 w l lw 2 title "SC", \
     "eam_exafs.chir" u 1:2 w l lw 2 title "EAM", \
     "meam_exafs.chir" u 1:2 w l lw 2 title "MEAM"
