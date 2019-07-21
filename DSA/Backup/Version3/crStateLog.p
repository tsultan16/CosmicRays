#CR State Plot
#plot from file: output.txt

set terminal gif animate delay 20 size 1280, 680
set output 'CRState.gif'

tSteps=150
nx=500
np=100



#set key font ",10"

set autoscale
do for [i=0:tSteps-1] {
  set title "Time Step = #".(i+1)
  set multiplot layout 1,3

  set logscale x 10
  #set logscale y 10

  set xlabel "p"
  set ylabel "log[p^4 f(p)]"
  #set xrange[1e-2:1e3]
  #set xrange [1e-2:1e2]
#  plot "fp.txt" every ::i*np+1::np+(i*np)-1 using 1:2 with linespoint pointtype 7 lc rgb "red" title "Numerical" ,\
#  "fp.txt" every ::i*np+1::np+(i*np)-1 using 1:3 with points pointtype 4 lc rgb "blue" title "Exact"
  plot "fp.txt" every ::i*np+1::np+(i*np)-1 using 1:2 with linespoint pointtype 6 lc rgb "red" notitle
  unset title
  unset xlabel
  unset ylabel
  #unset xrange
  #unset yrange
  
  unset logscale x
  #unset logscale y

  set xlabel "x"
  set ylabel "f(x)"
  #set xrange [0:1]
  #set yrange[0:1.] 
  plot "fx.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:2 with linespoint pointtype 7 lc rgb "red" notitle
  unset title
  unset xlabel
  unset ylabel
  #unset xrange
  #unset yrange 

  set xlabel "x"
  set ylabel "Fluid Density"
  plot "output_fluid.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:2 with linespoint pointtype 7 lc rgb "blue" notitle 
  unset title
  unset xlabel
  unset ylabel
  unset yrange

  #set xlabel "x"
  #set ylabel "Fluid Velocity"
  #plot "output_fluid.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:3 with linespoint pointtype 7 lc rgb "blue" notitle 
  #unset title
  #unset xlabel
  #unset ylabel
  #unset yrange


  unset multiplot
  
}

unset output