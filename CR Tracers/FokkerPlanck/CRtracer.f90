!CosmicRay Fokker-Plack solver in Lagrangian frame
!using passive tracer particles. (3 phase space co-ordinates: x,y,|p|)

program CRtracers
use constants_mod
use readTracer_module
use FokkerPlanck_mod

implicit none

integer ::i,j,k,tempcell(1)

real*8::dens,vel,pres,dudxtemp,tstep(nt)

real::start,finish
call cpu_time(start)

!Cosmic Ray initialization
call CRinit()

t=0.

open(unit=10,file='Input/dt.txt')
do i=1,nt
  read(10,fmt=*) tstep(i)
end do

!Compute Solution
do i=1,nt
  print*,'TIME STEP=',i

  !read time-step size from file
  dt=tstep(i)
  print*,'dt=',dt

  !set time-step size
  t=t+dt

  !read tracer data for current time-step 
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
  end if
  filename=trim('Input/tracer_t=')//trim(uniti)//trim('.txt')
  call read_data()
  
  
  print*,'Updating CR...'
  !Evolve Cosmic Ray Distribution Function 
  call CRevolve()  
  print*,'CR Update Done.'
  

end do 


print*,'Done.'

close(unit=11)
close(unit=20)

call cpu_time(finish)

print*,'Time Elapsed (seconds)=',finish-start

contains



end program CRTracers
