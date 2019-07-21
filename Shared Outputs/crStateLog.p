#CR State Plot
#plot from file: output.txt

set terminal gif animate delay 20 size 1280, 680
set output 'CRState.gif'

tSteps=100
nx=500
np1=50
np2=1000


#set key font ",10"

set autoscale
do for [i=0:tSteps-1] {
  set title "Time Step = #".(i+1)

  set logscale x 10
  set logscale y 10

  set xlabel "p"
  set ylabel "p^4 f(p)"
  #set xrange[1e-3:1e6]
  set yrange[1e-8:1e-4]

  plot "fp1.txt" every ::i*np1+1::np1+(i*np1)-1 using 1:2 with linespoint pointtype 6 lc rgb "red" title "CGMV" ,\
  "fp2.txt" every ::i*np2+1::np2+(i*np2)-1 using 1:2 with linespoint pointtype 6 lc rgb "blue" title "FD"
  unset title
  unset xlabel
  unset ylabel
  #unset xrange
  unset yrange
  
  unset logscale x
  unset logscale y


  
}

unset output
