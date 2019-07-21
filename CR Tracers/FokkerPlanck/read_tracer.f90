!Descrition: Contains subroutine that reads and stores data from "tvd.out" output file

module readTracer_module
use constants_mod

implicit none


integer,parameter ::lun=10
integer ::res
real*8 ::tstart, tstop, tdump
integer ::tSteps

character(len=6)::uniti
character(len=40)::filename


!file/io variables packed into a derived data type
type iovar
  !integer ::N
  !character(len=80) ::filename
  character(len=80) ::cbuffer
  integer ::cpos,strlen
  real*8,pointer::vel(:,:),rho(:),divU(:),mag(:,:)
  integer,pointer::id(:),cell(:,:)
end type iovar


!declare variable with the derived data type
type(iovar) file_dat

!subroutine for reading file contents
contains

subroutine read_data()

integer ::i

!open file (the variable following iostat= is set to a 
!non-zero integer if error is detected)
open(unit=lun,file=filename,form='formatted',iostat=res)
if(res/=0) then   
 close(unit=lun)
 print *,'Error opening file, status: ',res
 stop
end if
print*,'File open successful.'

!allocate memory for data arrays accordingly
allocate(file_dat%id(N)) !tracer id
allocate(file_dat%cell(N,3)) !cell-coordinate
allocate(file_dat%vel(N,3)) !fluid velocity
allocate(file_dat%rho(N))  !fluid density
allocate(file_dat%divU(N)) !fluid velocity divergence
!allocate(file_dat%mag(N,3)) !magneitc field

!start reading data from file and storing in arrays
do i=1,N
  read(unit=lun,fmt=*) file_dat%id(i) &
  ,file_dat%cell(i,1),file_dat%cell(i,2) &
  ,file_dat%vel(i,1),file_dat%vel(i,2) &
  ,file_dat%rho(i),file_dat%divU(i)
end do

!file contents have been read. Now, close file
close(unit=lun)
	
end subroutine read_data

end module readTracer_module
