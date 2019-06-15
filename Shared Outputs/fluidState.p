#Fluid State Plot
#plot from file: output.txt

set terminal gif animate delay 10 size 1280, 680
set output 'convergence_fluid.gif'


tSteps=100
nx1=500
nx2=nx1
nx3=nx2
nx4=nx3


set key font ",10"

set autoscale
do for [i=0:tSteps-1] {
  set title "Time Step = #".(i+1)
  set multiplot layout 3,1
  #plot 1: fluid density
  #set label 1 'Fluid Density' at graph 0.92,0.9 font ',8'
  #set xlabel "x"
  set ylabel "Density"
  plot "output_fluid1.txt" every ::i*nx1+1::nx1+(i*nx1)-1 using 1:2 with lines lc rgb "red" title "CGMV" ,\
 "output_fluid2.txt" every ::i*nx2+1::nx2+(i*nx2)-1 using 1:2 with lines lc rgb "blue" title "FD" 
 
  unset title
  unset xlabel
  unset ylabel

  #plot 3: Fluid Pressure
  #set label 1 'Fluid Pressure' at graph 0.92,0.9 font ',8'
  set xlabel "x"
  set ylabel "Pgas"
   plot "output_fluid1.txt" every ::i*nx1+1::nx1+(i*nx1)-1 using 1:4 with lines lc rgb "red" title "CGMV" ,\
 "output_fluid2.txt" every ::i*nx2+1::nx2+(i*nx2)-1 using 1:4 with lines lc rgb "blue" title "FD" 
 
  
  unset xlabel
  unset ylabel

#plot 3: CR Pressure
  #set label 1 'CR Pressure' at graph 0.92,0.9 font ',8'
  set xlabel "x"
  set ylabel "Pc"
  plot "output_fluid1.txt" every ::i*nx1+1::nx1+(i*nx1)-1 using 1:5 with lines lc rgb "red" title "CGMV" ,\
 "output_fluid2.txt" every ::i*nx2+1::nx2+(i*nx2)-1 using 1:5 with lines lc rgb "blue" title "FD" 
 
  
  unset xlabel
  unset ylabel


  unset multiplot

}

unset output
