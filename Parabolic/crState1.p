#CR State Plot
#plot from file: output.txt

set terminal gif animate delay 5 size 1280, 680
set output 'CRState.gif'

TOP=0.90
DY = 0.23
tSteps=200
np=100
nx=100



#set key font ",10"

set autoscale
do for [i=0:tSteps-1] {
  set title "Time Step = #".(i+1)
  set multiplot layout 2,1

  #set logscale x 10
  #set logscale y 10
  set yrange[0:26.5]

  set xlabel "p"
  set ylabel "f(p)"
  plot "fp.txt" every ::i*np+1::np+(i*np)-1 using 1:2 with linespoint pointtype 7 lc rgb "red" title "implicit" ,\
  "fp2.txt" every ::i*np+1::np+(i*np)-1 using 1:2 with linespoint pointtype 6 lc rgb "blue" title "semi-implicit" ,\
  "fp3.txt" every ::i*np+1::np+(i*np)-1 using 1:2 with linespoint pointtype 5 lc rgb "green" title "explicit"
  #set xrange [1e-3:1e3]
  unset title
  unset xlabel
  unset ylabel
  #unset xrange

  #unset logscale x
  #unset logscale y

set xlabel "x"
  set ylabel "f(x)"
  plot "fx.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:2 with linespoint pointtype 7 lc rgb "red" notitle
  set xrange [0.:1.] 
  unset title
  unset xlabel
  unset ylabel
  unset xrange 

  unset multiplot
  
}

unset output
