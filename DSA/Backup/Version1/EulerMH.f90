!Roe upwind method (1st order) for Eulers Equation in 1D
!with Cosmic Ray Fokker-Planck Solver

program EulerMH
use constants_mod
use RoeSolver_entropyfixed_mod
use FokkerPlanck_mod

implicit none

integer ::i,j,k
character(len=6)::uniti
character(len=20) :: filename

real*8::dens,vel,pres,shockPos,temp


real::start,finish
call cpu_time(start)

!Cosmic Raty and Fluid initialization
call init()
call CRinit()

open(unit=123,file='dudx.txt')


!Compute Solution
do i=1,nt
  print*,'TIME STEP=',i

  !compute numerical fluxes
  print*,'Updating Fluid...'
  call computeFlux()


  !set time-step size
  dt=dx*cour/smax

  print*,'dt=',dt
  !write(12,*) dt
  do j=1,nx   
    !update cell averages of conerved variables
    do k=1,3    
     u2(j,k)=u1(j,k)+(dt/dx)*(flux(j-1,k)-flux(j,k))
    end do


    x=xmin+(j-0.5)*dx 
    dens=u2(j,1)
    vel=u2(j,2)/u2(j,1)
    pres=(gam-1.)*(u2(j,3)-0.5*u2(j,2)*u2(j,2)/u2(j,1))
    if(mod(i,tSkip)==0)then
      write(11,*) x,dens,vel,pres
    end if
  end do  

  print*,'Fluid Update Done.'

  do j=1,nx
   dudx(j)=(u3(j,2)/u3(j,1)-u3(j-1,2)/u3(j-1,1))/dx
   if(mod(i,tSkip)==0)then
    write(123,*) xmin+(j-0.5)*dx,dudx(j)
   end if 
  end do

  print*,'Updating CR...'
  !Evolve Cosmic Ray Distribution Function  
  call CRevolve()  


  !Basic Shock Tracker
  temp=0.
  shockCell=0
  do j=1,nx
   !print*,'j,dudx=',j,dudx(j)
   shockPos=max(temp,abs(dudx(j)))
   if(shockPos .ne. temp .and. j>shockCell) shockCell=j
  end do
  print*,'Shock located near  cell#',shockCell 
  
  !print*,'CHECKPOINT1'
  if(mod(i,tSkip)==0)then
   do j=0,nx-1
    write(13,*) xmin+j*dx,f(50,j)!fx(j)
   end do
   do j=0,np-1
    write(10,*) px(j),f(j,5)!fp(j)
   end do
  end if
  !print*,'CHECKPOINT2' 

  print*,'CR Update Done.'

  !Enforce fluid boundary conditions
  call fluidBound()	

  u1=u2


end do 


print*,'Done.'
close(unit=11)
close(unit=10)
close(unit=13)
close(unit=123)

call cpu_time(finish)

print*,'Time Elapsed (seconds)=',finish-start

contains

end program EulerMH
