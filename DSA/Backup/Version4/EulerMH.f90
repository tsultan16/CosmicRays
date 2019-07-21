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

  !compute numerical fluxes
  print*,'Updating Fluid...'
  call computeFlux()
  
  !set time-step size
  dt=cour*dx/smax
  t=t+dt

  print*,'dt=',dt
  !write(12,*) dt
  do j=0,nx-1   
    !update cell averages of conerved variables
    do k=1,3    
     u2(j,k)=u1(j,k)+(dt/dx)*(flux(j-1,k)-flux(j,k))
    end do
 
    !----------------------------------------------------------------------------
    !Include CR pressure back reaction effect on fluid momentum 
    u2(j,2)=u2(j,2)-(dt/dx)*0.5*(Pc(j+1)-Pc(j-1))
    !Include CR pressure back reaction effect on fluid energy 
    u2(j,3)=u2(j,3)-(dt/dx)*0.5*(Pc(j+1)-Pc(j-1))*(u1(j,2)/u1(j,1))
    u2(j,3)=u2(j,3)-dt*fluidS(j)*dx
    !----------------------------------------------------------------------------

    x=xmin+(j-0.5)*dx 
    dens=u2(j,1)
    vel=u2(j,2)/u2(j,1)
    pres=(gam-1.)*(u2(j,3)-0.5*u2(j,2)*u2(j,2)/u2(j,1))
    if(mod(i,tSkip)==0)then
      write(11,*) x,dens,vel,pres,Pc(j)
    end if
  end do  


  !Enforce fluid boundary conditions
  call fluidBound()

  print*,'Fluid Update Done.'

  do j=0,nx-1
   !dudx(j)=(u3(j,2)/u3(j,1)-u3(j-1,2)/u3(j-1,1))/dx
   dudx(j)=(u1(j,2)/u1(j,1)-u1(j-1,2)/u1(j-1,1))/dx
   if(mod(i,tSkip)==0)then
    write(124,*) xmin+(j-0.5)*dx,dudx(j)
   end if 
  end do


  !Basic Shock Tracker (locates cell with largest velocity gradient)
  shockCell=minloc(dudx)
  print*,'Shock located near  cell#',shockCell(1)
  !dudxtemp=dudx(shockCell(1)-1)
  !dudx=0.
  !dudx(shockCell(1)-1)=dudxtemp

  do j=0,nx-1
   if(mod(i,tSkip)==0)then
    write(123,*) xmin+(j-0.5)*dx,dudx(j)
   end if 
  end do

  print*,'Updating CR...'
  !Evolve Cosmic Ray Distribution Function 
  call CRevolve()  
  
  !print*,'CHECKPOINT1'
  if(mod(i,tSkip)==0)then
   do j=0,nx-1
    write(13,*) xmin+j*dx,f(25,j)!f(1,j)!fx(j)
   end do
   do j=0,np-1
    write(10,*) px(j),(px(j)**4.)*f(j,shockCell(1)-2),fexact(j)!fp(j)
   end do
  end if
  !print*,'CHECKPOINT2' 

  print*,'CR Update Done.'

  u1=u2


end do 


print*,'Done.'
close(unit=11)
close(unit=10)
close(unit=13)
close(unit=123)
close(unit=124)

call cpu_time(finish)

print*,'Time Elapsed (seconds)=',finish-start

contains

end program EulerMH
