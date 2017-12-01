# Gnuplot script file for plotting data in file
set multiplot layout 2,2 rowsfirst
set macros
POS = "at graph 0.92,0.9 font ',8'"
NOXTICS = "set xtics ('' 0,'' 1000,'' 2000,'' 3000,'' 4000); \
          unset xlabel"
XTICS = "set xtics ('0' 0,'1000' 1000,'2000' 2000,'3000' 3000,'4000' 4000);\
          set xlabel 'Steps'"
NOYTICS = "set format y ''; unset ylabel"
YTICS = "set format y '%.1f'; set ylabel 'Energy(eV)'"

# Margins for each row resp. column
TMARGIN = "set tmargin at screen 0.95; set bmargin at screen 0.65"
MMARGIN = "set tmargin at screen 0.65; set bmargin at screen 0.35"
BMARGIN = "set tmargin at screen 0.35; set bmargin at screen 0.05"
LMARGIN = "set lmargin at screen 0.15; set rmargin at screen 0.55"
RMARGIN = "set lmargin at screen 0.55; set rmargin at screen 0.95"

set label 1 'a' @POS
@NOXTICS; @YTICS; @TMARGIN; @LMARGIN
plot[0:5000][-190:-186] "0_02.dat" using 2:7 with lines

@NOXTICS; @NOYTICS; @TMARGIN; @RMARGIN
plot[0:5000][-190:-186] "0_04.dat" using 2:7 with lines

@NOXTICS; @YTICS; @MMARGIN; @LMARGIN
plot[0:5000][-190:-186] "0_08.dat" using 2:7 with lines

@NOXTICS; @NOYTICS; @MMARGIN; @RMARGIN
plot[0:5000][-190:-186] "0_16.dat" using 2:7 with lines

@XTICS; @YTICS; @BMARGIN; @LMARGIN
plot[0:5000][-190:-186] "0_32.dat" using 2:7 with lines

@XTICS; @NOYTICS; @BMARGIN; @RMARGIN
plot[0:5000][-190:-186] "0_50.dat" using 2:7 with lines
