#CR State Plot
#plot from file: output.txt

set terminal gif animate delay 10 size 1280, 680
set output 'CRState.gif'

tSteps=100
np=100



#set key font ",10"

set autoscale
do for [i=0:tSteps-1] {
  set title "Time Step = #".(i+1)

  set logscale x 10
  set logscale y 10

  set xlabel "p"
  set ylabel "p^4 f(p)"
  #set xrange[1e-3:1e6]
  set yrange[1e-14:1e-1]

  plot "fp1.txt" every ::i*np+1::np+(i*np)-1 using 1:2 with linespoint pointtype 6 lc rgb "red" title "dt2=1.5*dt" ,\
 "fp2.txt" every ::i*np+1::np+(i*np)-1 using 1:2 with linespoint pointtype 6 lc rgb "green" title "dt2=2.0*dt" ,\
 "fp3.txt" every ::i*np+1::np+(i*np)-1 using 1:2 with linespoint pointtype 6 lc rgb "blue" title "dt2=2.5*dt" ,\
 "fp4.txt" every ::i*np+1::np+(i*np)-1 using 1:2 with linespoint pointtype 6 lc rgb "purple" title "dt2=3.0*dt"
  
  unset title
  unset xlabel
  unset ylabel
  #unset xrange
  unset yrange
  
  unset logscale x
  unset logscale y



  unset multiplot
  
}

unset output
