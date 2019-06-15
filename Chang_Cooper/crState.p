#CR State Plot
#plot from file: output.txt

set terminal gif animate delay 10 size 1280, 680
set output 'CRState.gif'

TOP=0.90
DY = 0.23
tSteps=400
nx=100

set logscale x 10
set logscale y 10

#set key font ",10"

set autoscale
do for [i=0:tSteps-1] {
  set title "Time Step = #".(i+1)
  set xlabel "p"
  set ylabel "f(p)"
  plot "fp.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:2 with linespoint pointtype 7 lc rgb "red" title "Implicit" ,\
  "fp2.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:2 with linespoint pointtype 7 lc rgb "blue" title "Semi Implicit"
  set xrange [0.0001:10000]
  unset title
  unset xlabel
  unset ylabel
  unset yrange
  
}

unset output
