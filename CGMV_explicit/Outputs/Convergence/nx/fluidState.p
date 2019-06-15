#Fluid State Plot
#plot from file: output.txt

set terminal gif animate delay 30 size 1280, 680
set output 'convergence_fluid.gif'


tSteps=100
nx1=500
nx2=1000

tSkip=2

set key font ",10"

set autoscale
do for [i=0:tSteps-1] {
  j=i*tSkip
  set title "Time Step = #".(i+1)
  set multiplot layout 3,1
  
  set xlabel "x"
  set ylabel "Density"
  plot "output_fluid2.txt" every ::j*nx2+1::nx2+(j*nx2)-1 using 1:2 with lines lc rgb "red" title "nx=1000" ,\
  "output_fluid4.txt" every ::j*nx2+1::nx2+(j*nx2)-1 using 1:2 with lines lc rgb "green" title "nx=1000"  
  unset title
  unset xlabel
  unset ylabel

  #plot 3: Fluid Pressure

  set xlabel "x"
  set ylabel "Pgas"
  plot "output_fluid2.txt" every ::j*nx2+1::nx2+(j*nx2)-1 using 1:4 with lines lc rgb "red" title "nx=1000" ,\
  "output_fluid4.txt" every ::j*nx2+1::nx2+(j*nx2)-1 using 1:4 with lines lc rgb "green" title "nx=1000"
  unset xlabel
  unset ylabel

#plot 3: CR Pressure

  set xlabel "x"
  set ylabel "Pc"
    plot "output_fluid2.txt" every ::j*nx2+1::nx2+(j*nx2)-1 using 1:5 with lines lc rgb "red" title "nx=1000" ,\
  "output_fluid4.txt" every ::j*nx2+1::nx2+(j*nx2)-1 using 1:5 with lines lc rgb "green" title "nx=1000"
  unset xlabel
  unset ylabel


  unset multiplot

}

unset output
