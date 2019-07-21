!Van Leer MUSCL Scheme (2nd order) for Ieal MHD in 2D
!With CT for maintaining div(B)=0 on the grid

program MHD_MH
use constants_mod
use RoeSolver2D_mod
use velocity_tracer_mod_2D
use mass_tracer_mod_2D

implicit none

integer::i,j,k,l,alt_flag
character(len=6)::uniti
real::start,finish
real*8::t,yx(nx),yx2(nx),xt0,yt0,yt2_0!,xt(0:nx),yt(0:ny),xt2(0:nx),yt2(0:ny)
real*8::rad


t=0._8


open(unit=99,file='Output/energy.txt')
open(unit=88,file='Output/dt.txt')



call cpu_time(start)

!fluid initialization
call init()
!tracer initialization
if(tracerType==1)then
 call MassTracerInit()
else if(tracerType==2)then
 call velocityTracerInit()
end if

alt_flag=1

!Compute Solution 
do i=1,nt
  call timeStep() 
  print*,'timestep,dt',i,dt
  t=t+dt
  write(88,*) dt

  print*,'Updating fluid.'
  if(alt_flag==0)then

    call x_sweep()
    call y_sweep()
    alt_flag=1
  else if(alt_flag==1)then
    call y_sweep()
    call x_sweep()
    alt_flag=0
  end if

  !First Order Correction for uniform gravitaional field (source term)
  if(gravity==1)then
  do j=1,nx 
   x=xmin+(j-0.5)*dx
   do k=1,ny
     y=ymin+(k-0.5)*dy
     rad=sqrt((x-0.5)**2+(y-0.5)**2)
     gx=-0.5*((x-0.5)/rad)
     gy=-0.5*((y-0.5)/rad)

     u2(j,k,8)=u2(j,k,8)+dt*(u2(j,k,2)*gx+u2(j,k,3)*gy)
     !u2(j,k,8)=u2(j,k,8)-0.5*(u2(j,k,2)**2.+u2(j,k,3)**2.+u2(j,k,4)**2.)/u2(j,k,1) 
     u2(j,k,2)=u2(j,k,2)+dt*u2(j,k,1)*gx !x-momentum
     u2(j,k,3)=u2(j,k,3)+dt*u2(j,k,1)*gy !y-momentum
     !u2(j,k,8)=u2(j,k,8)+0.5*(u2(j,k,2)**2.+u2(j,k,3)**2.+u2(j,k,4)**2.)/u2(j,k,1) !total fluid energy
   end do
  end do
  call bound()
  u1=u2
  end if

   
  !Update Cell interface Magnetic Field (CT)
  call magUpdate()


  call horcut()
  !call vercut()
  !call diagcut()
   
  !Compute Velocity Divergence
  do j=1,nx 
   do k=1,ny
   divU(j,k)=(u2(j,k,2)/u2(j,k,1)-u2(j-1,k,2)/u2(j-1,k,1))/dx &
             +(u2(j,k,3)/u2(j,k,1)-u2(j,k-1,3)/u2(j,k-1,1))/dy
   end do
  end do

  print*,'Fluid update complete.'  
  !call computeEnergy()

  if(tracerType==2)then
    call velTracerAdvect(dt)
  end if


  if(mod(i,tSkip)==0)then

  if(i<10)then
    write(uniti,'(I1.1)') i
  else if(i>=10 .and. i<100)then
    write(uniti,'(I2.2)') i
  else if(i>=100 .and. i<1000)then
    write(uniti,'(I3.3)') i
  else if(i>=1000 .and. i<10000)then
    write(uniti,'(I4.3)') i
  else if(i>=10000 .and. i<100000)then
    write(uniti,'(I5.3)') i
  else if(i>=100000 .and. i<1000000)then
    write(uniti,'(I6.3)') i
  end if
  filename1=trim('Output/t=')//trim(uniti)//trim('.txt')
  filename2=trim('Output/tracer_t=')//trim(uniti)//trim('.txt')

  open(unit=i+1000,file=filename1)
  open(unit=i+20000,file=filename2)
  

  if(tracerType==1)then
   call fileOutput(i+1000,N_cell)
   call massTracerOutput(i+20000)  
  else if(tracerType==2)then
   call fileOutput(i+1000,cell_count)
   call velocityTracerOutput(i+20000) 
  end if
 
  close(unit=i+1000)
  close(unit=i+20000)
 
  end if
  
  
  print*,'% Complete=',(real(i)/real(nt))*100.


end do 

close(unit=12)
!close(unit=13)
close(unit=14)
close(unit=15)
close(unit=16)
close(unit=17)

close(unit=20)
close(unit=88)
close(unit=99)

print*,'Done.'
call cpu_time(finish)
print*,'Time Elapsed (seconds)=',finish-start


contains

!-----------------------------------------------
!x-sweeps
!----------------------------------------------- 
subroutine x_sweep() 
  
do k=0,ny+1
 !compute numerical fluxes
 if(fluxType==1)then    
   call computeFluxRoe(nx,k,1)
 else if(fluxType==2)then
   call computeFluxHLL(nx,k,1)
 else if(fluxType==3)then
   call computeFluxHLLI(nx,k,1)
 end if
 do j=1,nx   
   !update cell averages of conserved variables
   u2(j,k,1)=u1(j,k,1)+(dt/dx)*(flux(j-1,1)-flux(j,1)) 
   u2(j,k,2)=u1(j,k,2)+(dt/dx)*(flux(j-1,2)-flux(j,2))
   u2(j,k,3)=u1(j,k,3)+(dt/dx)*(flux(j-1,3)-flux(j,3))
   u2(j,k,4)=u1(j,k,4)+(dt/dx)*(flux(j-1,4)-flux(j,4))
   !u2(j,k,5)=u1(j,k,5)
   !u2(j,k,6)=u1(j,k,6)+(dt/dx)*(flux(j-1,5)-flux(j,5))
   u2(j,k,7)=u1(j,k,7)+(dt/dx)*(flux(j-1,6)-flux(j,6))
   u2(j,k,8)=u1(j,k,8)+(dt/dx)*(flux(j-1,7)-flux(j,7)) 

   !Mass tracer advection
   if(tracerType==1)then 
    call massTracerAdvect(j,k,flux(j-1,1),flux(j,1),u1(j,k,1),dt,dx,1) 
   end if

 end do
end do

call protection()
call bound()	
u1=u2

end subroutine x_sweep 

!-----------------------------------------------
!y-sweeps
!-----------------------------------------------
subroutine y_sweep()

do j=0,nx+1
 !compute numerical fluxes
 if(fluxType==1)then    
   call computeFluxRoe(ny,j,2)
 else if(fluxType==2)then    
   call computeFluxHLL(ny,j,2)
 else if(fluxType==3)then    
   call computeFluxHLLI(ny,j,2)
 end if 
 do k=1,ny   
   !update cell averages of conserved variables    
   u2(j,k,1)=u1(j,k,1)+(dt/dy)*(flux(k-1,1)-flux(k,1)) 
   u2(j,k,2)=u1(j,k,2)+(dt/dy)*(flux(k-1,3)-flux(k,3))
   u2(j,k,3)=u1(j,k,3)+(dt/dy)*(flux(k-1,2)-flux(k,2))
   u2(j,k,4)=u1(j,k,4)+(dt/dy)*(flux(k-1,4)-flux(k,4))
   !u2(j,k,5)=u1(j,k,5)+(dt/dy)*(flux(k-1,5)-flux(k,5))
   !u2(j,k,6)=u1(j,k,6)
   u2(j,k,7)=u1(j,k,7)+(dt/dy)*(flux(k-1,6)-flux(k,6))
   u2(j,k,8)=u1(j,k,8)+(dt/dy)*(flux(k-1,7)-flux(k,7))  

   !Mass tracer advection 
   if(tracerType==1)then
    call massTracerAdvect(j,k,flux(k-1,1),flux(k,1),u1(j,k,1),dt,dy,2) 
   end if
    
 end do
end do

call protection()
call bound()	
u1=u2

end subroutine y_sweep

subroutine magUpdate()

!------------------------------------------
!Cell Corner Electric Field: Ez[i,j]
!------------------------------------------
do j=0,nx
!$omp parallel do shared(Ez,fx,fy)
  do k=0,ny
   Ez(j,k)=0.5*(fy(j,k)+fy(j+1,k)-fx(j,k)-fx(j,k+1))
  end do
 !$omp end parallel do  
end do

!------------------------------------------------------------------------
!Update Cell interface magnetic field normal components: bx[i,j], by[i,j]
!-----------------------------------------------------------------------
do k=1,ny
  !$omp parallel do shared(Ez,bxInt)
  do j=0,nx
   bxInt(j,k)=bxInt(j,k)-(dt/dy)*(Ez(j,k)-Ez(j,k-1))
  end do
  !$omp end parallel do   
end do

do j=1,nx
  !$omp parallel do shared(Ez,byInt)
  do k=0,ny
   byInt(j,k)=byInt(j,k)+(dt/dx)*(Ez(j,k)-Ez(j-1,k))
  end do
   !$omp end parallel do    
end do


!-----------------------------------------------------------------
!Update Cell Center Magnetic Field Values (linear interpolation 
!of cell inteferface values)
!-----------------------------------------------------------------
do j=1,nx
  !$omp parallel do shared(u2,bxInt,byInt)
  do k=1,ny
    u2(j,k,5)=0.5*(bxInt(j-1,k)+bxInt(j,k))
    u2(j,k,6)=0.5*(byInt(j,k-1)+byInt(j,k))
  end do
  !$omp end parallel do
end do

call bound()
u1=u2


end subroutine magUpdate


subroutine computeEnergy()
real*8::eint,ekin,emag

eint=0.
ekin=0.
emag=0.

do j=1,nx 
 do k=1,ny
   ekin=ekin+dx*dy*0.5*(u2(j,k,2)**2.+u2(j,k,3)**2.+u2(j,k,4)**2.)/u2(j,k,1)
   emag=emag+dx*dy*0.5*(u2(j,k,5)**2.+u2(j,k,6)**2.+u2(j,k,7)**2.)
   eint=eint+dx*dy*( u2(j,k,8)&
   -0.5*(u2(j,k,2)**2.+u2(j,k,3)**2.+u2(j,k,4)**2.)/u2(j,k,1) &
   -0.5*(u2(j,k,5)**2.+u2(j,k,6)**2.+u2(j,k,7)**2.) )
 end do
end do

write(99,*) t,ekin,eint,emag 

end subroutine computeEnergy


end program MHD_MH
