#Fluid State Plot
#plot from file: output.txt

set terminal gif animate delay 20 size 1280, 680
set output 'FluidState.gif'


tSteps=200
nx=1000

set key font ",10"

set autoscale
do for [i=0:tSteps-1] {
  set title "Time Step = #".(i+1)
  set multiplot layout 2,2
  #plot 1: fluid density
  #set label 1 'Fluid Density' at graph 0.92,0.9 font ',8'
  #set xlabel "x"
  set ylabel "Density"
  plot "output_fluid.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:2 with lines lc rgb "blue" notitle 
  unset title
  unset xlabel
  unset ylabel

  #plot 2: Fluid Velocity
  set xlabel "x"
  set ylabel "Fluid Velocity"
  plot "output_fluid.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:3 with lines lc rgb "blue" notitle 
  unset xlabel
  unset ylabel
 
  #plot 3: Fluid Pressure
  #set label 1 'Fluid Pressure' at graph 0.92,0.9 font ',8'
  set xlabel "x"
  set ylabel "Pgas"
  plot "output_fluid.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:4 with lines lc rgb "blue" notitle  
  unset xlabel
  unset ylabel

#plot 3: CR Pressure
  #set label 1 'CR Pressure' at graph 0.92,0.9 font ',8'
  set xlabel "x"
  set ylabel "Pc"
  plot "output_fluid.txt" every ::i*nx+1::nx+(i*nx)-1 using 1:5 with lines lc rgb "blue" notitle  
  unset xlabel
  unset ylabel


  unset multiplot

}

unset output
