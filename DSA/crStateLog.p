#CR State Plot
#plot from file: output.txt

set terminal gif animate delay 100 size 1280, 680
set output 'CRState.gif'

tSteps=100
nx=1000
np=500
tSkip=1


#set key font ",10"

set autoscale
do for [i=0:tSteps-1] {
  j=i*tSkip
  set title "Time Step = #".(i+1)
  set multiplot layout 1,2

  set logscale x 10
  set logscale y 10

  set xlabel "p"
  set ylabel "log[p^4 f(p)]"
  #set xrange[1e-2:1e20]
  set yrange[1e-12:1]

  plot "fp.txt" every ::j*np+1::np+(j*np)-1 using 1:2 with linespoint pointtype 6 lc rgb "red" notitle
  unset title
  unset xlabel
  unset ylabel
  #unset xrange
  unset yrange
  
  unset logscale x
  unset logscale y

  set xlabel "x"
  set ylabel "f(x)"
  #set xrange [0:1]
  #set yrange[0:1.] 
  plot "fx.txt" every ::j*nx+1::nx+(j*nx)-1 using 1:2 with linespoint pointtype 7 lc rgb "red" notitle
  unset title
  unset xlabel
  unset ylabel
  #unset xrange
  unset yrange 



  unset multiplot
  
}

unset output
