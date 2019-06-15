!Roe upwind method (1st order) for Eulers Equation in 1D
!with Cosmic Ray Fokker-Planck Solver

program EulerMH
use constants_mod
use RoeSolver_entropyfixed_mod
use FokkerPlanck_mod

implicit none

integer ::i,j,k,tempcell(1)
character(len=6)::uniti
character(len=20) :: filename

real*8::dens,vel,pres,dudxtemp

real::start,finish
call cpu_time(start)

!Cosmic Raty and Fluid initialization
call init()
call CRinit()

open(unit=123,file='dudx.txt')
open(unit=124,file='dudx_original.txt')

t=0.

!Compute Solution
do i=1,nt
  print*,'TIME STEP=',i

   !set time-step size
  dt=0.01
  t=t+dt

  
  !Enforce fluid boundary conditions
  call fluidBound()

  print*,'Fluid Update Done.'

  print*,'Updating CR...'
  !Evolve Cosmic Ray Distribution Function 
  call CRevolve()    

  if(mod(i,tSkip)==0)then
   do j=0,nx-1
    write(13,*) xmin+j*dx,f(25,j)
   end do
   do j=0,np-1
    write(10,*) px(j),f(j,1)
   end do
  end if


  print*,'CR Update Done.'
 
  u1=u2


end do 


print*,'Done.'
close(unit=11)
close(unit=10)
close(unit=13)
close(unit=123)
close(unit=124)

close(unit=110)
close(unit=111)


call cpu_time(finish)

print*,'Time Elapsed (seconds)=',finish-start

contains

end program EulerMH
