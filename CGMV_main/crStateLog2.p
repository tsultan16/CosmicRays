#CR State Plot

set terminal gif animate delay 10 size 1280, 680
set output 'CRState2.gif'

tSteps=300
np=20
np=np-1
nx=1000

tSkip=1

#set key font ",10"

set autoscale
do for [i=0:tSteps-1] {
  j=i*tSkip
  set title "Time Step = #".(i+1)
  set multiplot layout 2,3

  set logscale x 10
  #set logscale y 10

  set xlabel "p"
  set ylabel "n(p)"
  #set xrange[1e-3:1e6]
  #set yrange[0:1.3e6]

  plot "ng_p.txt" every ::j*np+1::np+(j*np)-1 using 2:3 with linespoint pointtype 5 lc rgb "red" notitle
  unset title
  unset xlabel
  unset ylabel
  #unset xrange
  #unset yrange

  set xlabel "p"
  set ylabel "g(p)"
  #set xrange[1e-3:1e6]
  #set yrange[1.e-6:5e-3]

  plot "ng_p.txt" every ::j*np+1::np+(j*np)-1 using 2:4 with linespoint pointtype 5 lc rgb "red" notitle
  unset title
  unset xlabel
  unset ylabel
  #unset xrange
  #unset yrange

  #unset logscale y

  set xlabel "p"
  set ylabel "q(p)"
  #set xrange[1e-3:1e6]
  #set yrange[2:6]
  plot "ng_p.txt" every ::j*np+1::np++(j*np)-1 using 2:5 with linespoint pointtype 5 lc rgb "red" notitle
  unset title
  unset xlabel
  unset ylabel
  #unset xrange
  unset yrange

  unset logscale x


  set xlabel "x"
  set ylabel "n(x)"
  plot "ng_x.txt" every ::j*nx+1::nx+(j*nx)-1 using 2:3 with lines lc rgb "red" notitle
  unset title
  unset xlabel
  unset ylabel


  set xlabel "x"
  set ylabel "g(x)"
  plot "ng_x.txt" every ::j*nx+1::nx+(j*nx)-1 using 2:4 with lines lc rgb "red" notitle
  unset title
  unset xlabel
  unset ylabel


  set xlabel "x"
  set ylabel "q(x)"
  #set yrange[2:6]
  plot "ng_x.txt" every ::j*nx+1::nx+(j*nx)-1 using 2:5 with linespoint pointtype 5 lc rgb "red" notitle
  unset title
  unset xlabel
  unset ylabel
  #unset yrange


  unset multiplot
  
}

unset output
