#CR State Plot
#plot from file: output.txt

set terminal gif animate delay 1 size 1280, 680
set output 'CRState.gif'

tSteps=500
np=500
tSkip=1


#set key font ",10"

print "Total time steps= ".(tSteps*tSkip)

set autoscale
do for [i=0:tSteps-1] {
  j=i*tSkip
  print "Time Steps Completed =".j

  set title "Time Step = #".(i+1)

  set logscale x 10
  set logscale y 10

  set xlabel "p"
  set ylabel "log[f(p)]"
  #set xrange[1e-2:1e20]
  set yrange[1e-12:1e5]

  plot "f(p)_1.txt" every ::j*np+1::np+(j*np)-1 using 1:2 with linespoint pointtype 6 lc rgb "red" notitle
  unset title
  unset xlabel
  unset ylabel
  #unset xrange
  unset yrange
  
  unset logscale x
  unset logscale y

  
}

unset output
