
for i in 0.0025 0.005 0.0075 0.01 0.0125 0.015
do
 grep 'free  energy   TOTEN' $i/OUTCAR|tail -1 >>y.dat
done
