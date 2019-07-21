module mass_tracer_mod_2D
use constants_mod
use RoeSolver2D_mod
implicit none



!---------------------------------------------------------------------------------
!Subroutine List:
!
!(1.) massTracerInit(u1,dx,dy),tracerDensity(xmin,ymin,dx,dy),
!(2.) cellTransfer(ix0,iy0,ix1,iy1,N_out),
!(3.) massTracerAdvect(ix,iy,fL,fR,rho,dt,dr,dir) 
!--------------------------------------------------------------------------------



!derived data types (i.e. struct ala c) for tracers
type mtr
  !real*8 ::mass,charge,current_pos(3),current_vel(3);
  real*8::velocity(3),density,pressure,divU,mag(3)
  integer ::cell_occupied(3) !cell occupied by the tracer
  integer ::tracer_id !assigned by order of "creation"
  
  integer::shockCell !set to 1 if the cell contains a shock
  real*8::shR !shock compression ratio
  real*8::shM !upstream mach number

  !integer ::tracer_role !0: passive, 1: active(useless for now)
  type(mtr),pointer ::next   !pointer for linked list capability
end type 
 
!cell_pointer points to head tracer (if any) in each cell 
type cell_ptr
  type(mtr), pointer ::p
end type cell_ptr

!tracer_pointer points to a tracer object
type tracer_ptr
  type(mtr), pointer ::p
end type tracer_ptr


real*8::tot_fluid_mass
real*8::m_cell !fluid mass in given cell

integer::N_cell(-1:nx+2,-1:ny+2) !number of tracers in a given cell
!temp variables
integer::N_cell_temp(-1:nx+2,-1:ny+2)
integer::temp_id

!create tracer_pointer and cell_pointers
type(tracer_ptr)::mTracer(N) !each array element is a pointer to a tracer

type(cell_ptr)::cell_tracer_head(0:nx+1,0:ny+1) !array will has nx elements (equal to number of cells), 
                                      !each element is a pointer to the "head" tracer 
  

type(mtr),pointer::tracer_temp, tracer_temp2, tracer_temp3 


contains

!--------------------------------------------------------------------
subroutine massTracerInit()
!--------------------------------------------------------------------
!Mass Tracer Initialization Routine

!Input Variables

!Local Variables
integer::ii,jj,kk,ixx,iyy
real*8::rand_num
integer::tr_counter,rand_cell(2)
real*8 ::tot_fluid_mass
real*8 ::m_cell !fluid mass in given cell

!open(unit=15,file='massTracerDensity.txt')
!open(unit=14,file='massTracerData.txt')


 
tr_counter=0
kk=1
N_cell=0
N_cell_temp=0
tot_fluid_mass=0._8

!allocate memory for first tracer
allocate(tracer_temp) 
tracer_temp%tracer_id =1
mTracer(1)%p=>tracer_temp

!----------------------------------------------------------------
!Match initial tracer distribution to initial fluid distribution
!----------------------------------------------------------------
if(init_option==1)then
 !compute total fluid mass
 do ii=1,nx
  do jj=1,ny
   tot_fluid_mass = tot_fluid_mass+u1(ii,jj,1)*dx*dy
  end do 
 end do

 do ii=1,nx
  do jj=1,ny
   !compute number of tracers to pe placed in jth cell
   m_cell = u1(ii,jj,1)*dx*dy !initial fluid mass in jth cell
   rand_num=rand(0)
   if(rand_num<0.5) then
     N_cell(ii,jj) = floor((m_cell/tot_fluid_mass)*N)
   else
     N_cell(ii,jj) = ceiling((m_cell/tot_fluid_mass)*N)
   end if 
   tr_counter = tr_counter+N_cell(ii,jj)
  end do
 end do
 
 !now distribute "leftover" tracers randomly across cells
 do while(tr_counter < N) 
   do ii=1,nx
    do jj=1,ny
     rand_cell(1)=nx*rand(0)
     rand_cell(2)=ny*rand(0)
     N_cell(rand_cell(1),rand_cell(2))=N_cell(rand_cell(1),rand_cell(2))+1 
     tr_counter = tr_counter+1
    end do
   end do
 end do
 !if excess tracers, subtract off the excess amount from cells selected randomly
 do while(tr_counter>N)
   rand_cell(1)=nx*rand(0)
   rand_cell(2)=ny*rand(0)
   if(N_cell(rand_cell(1),rand_cell(2))>1)then
    N_cell(rand_cell(1),rand_cell(2))=N_cell(rand_cell(1),rand_cell(2))-1
    tr_counter=tr_counter-1   
   end if
 end do
 

!----------------------------------------------------------------
!Uniformly distribute tracers in cells [minCell,maxCell]
!----------------------------------------------------------------
else if(init_option==2)then
 if((N<(maxCellx-minCellx+1)*(maxCelly-minCelly+1)))then
  print*,'Error: N too small for uniform distribution of tracers in [xmin,xmax]x[ymin,ymax]'
  stop
 end if
 do ii=minCellx,maxCellx
  do jj=minCelly,maxCelly
    !distribute uniformly across cells [minCellx,maxCellx]*[minCelly,maxCelly]
    !compute number of tracers to pe placed in ith cell
    rand_num=rand(0)
    if(rand_num<0.5) then
      N_cell(ii,jj) = floor(N*1./(1.*(maxCellx-minCellx+1)*(maxCelly-minCelly+1)))
    else
     N_cell(ii,jj) = ceiling(N*1./(1.*(maxCellx-minCellx+1)*(maxCelly-minCelly+1)))
    end if 
    tr_counter = tr_counter+N_cell(ii,jj)	  	  
  end do
 end do
 !now distribute "leftover" tracers uniformaly across cells
 do while(tr_counter < N) 
  do ii=minCellx,maxCellx
   do jj=minCelly,maxCelly
     N_cell(ii,jj)=N_cell(ii,jj)+1
     tr_counter = tr_counter+1
     if(tr_counter==N)then
      exit
     end if
   end do
  end do
 end do

 !if excess tracers, subtract off the excess amount from cells selected randomly
 do while(tr_counter>N)
   rand_cell(1)=minCellx+(maxCellx-minCellx+1)*rand(0)
   rand_cell(2)=minCelly+(maxCelly-minCelly+1)*rand(0)
   N_cell(rand_cell(1),rand_cell(2))=N_cell(rand_cell(1),rand_cell(2))-1
   tr_counter=tr_counter-1   
 end do

!-------------------------------------------------------------------------------------------
!Match initial tracer distribution to initial fluid distribution in cells [minCell,maxCell]
!-------------------------------------------------------------------------------------------
else if(init_option==3)then
 !compute total fluid mass
 do ii=minCellx,maxCellx
  do jj=minCelly,maxCelly
   tot_fluid_mass = tot_fluid_mass+u1(ii,jj,1)*dx*dy
  end do 
 end do

 do ii=minCellx,maxCellx
  do jj=minCelly,maxCelly
   !compute number of tracers to pe placed in jth cell
   m_cell = u1(ii,jj,1)*dx*dy !initial fluid mass in jth cell
   rand_num=rand(0)
   if(rand_num<0.5) then
     N_cell(ii,jj) = floor((m_cell/tot_fluid_mass)*N)
   else
     N_cell(ii,jj) = ceiling((m_cell/tot_fluid_mass)*N)
   end if 
   tr_counter = tr_counter+N_cell(ii,jj)
  end do
 end do
 
 !now distribute "leftover" tracers randomly across cells
 do while(tr_counter < N) 
   do ii=1,nx
    do jj=1,ny
     rand_cell(1)=nx*rand(0)
     rand_cell(2)=ny*rand(0)
     N_cell(rand_cell(1),rand_cell(2))=N_cell(rand_cell(1),rand_cell(2))+1 
     tr_counter = tr_counter+1
    end do
   end do
 end do
 !if excess tracers, subtract off the excess amount from cells selected randomly
 do while(tr_counter>N)
   rand_cell(1)=nx*rand(0)
   rand_cell(2)=ny*rand(0)
   if(N_cell(rand_cell(1),rand_cell(2))>1)then
    N_cell(rand_cell(1),rand_cell(2))=N_cell(rand_cell(1),rand_cell(2))-1
    tr_counter=tr_counter-1   
   end if
 end do
 



end if


print*,'tr_counter=',tr_counter

!place tracers in cells
do ixx=1,nx
 do iyy=1,ny
  if(N_cell(ixx,iyy)>0) then 
     !print*,'ixx,iyy,***dx,***dy,N_cell=',ixx,iyy,dx,dy,N_cell(ixx,iyy)
   !designate a head tracer for the (ix,iy)th cell, create pointer
   cell_tracer_head(ixx,iyy)%p=>mTracer(kk)%p !pointer to ith cell head tracer    
   mTracer(kk)%p%cell_occupied(1)=ixx
   mTracer(kk)%p%cell_occupied(2)=iyy
   !write tracer location to file
   !write(14,*) kk,mTracer(kk)%p%cell_occupied(1),mTracer(kk)%p%cell_occupied(2)
   !place tracers in (ix,iy)th cell, link head tracer to the others  
   if(N_cell(ixx,iyy)>1)then
     do jj=1,N_cell(ixx,iyy)-1	              
      !allocate memory for next tracer in the link. Tracer(k) is the "head tracer"	   
      allocate(mTracer(kk+jj-1)%p%next)
      mTracer(kk+jj)%p=>mTracer(kk+jj-1)%p%next
      !assign id
      mTracer(kk+jj)%p%tracer_id=kk+jj		
      !record cell location of this tracer
      mTracer(kk+jj)%p%cell_occupied(1)=ixx
      mTracer(kk+jj)%p%cell_occupied(2)=iyy
      !write tracer location to file**********
      !write(14,*) kk+jj,mTracer(kk+jj)%p%cell_occupied(1),mTracer(kk+jj)%p%cell_occupied(2) 
	
      if(jj==N_cell(ixx,iyy)-1) then	      
	!allocate memory for the head tracer of the next cell
        allocate(mTracer(kk+jj)%p%next)
	mTracer(kk+jj+1)%p=>mTracer(kk+jj)%p%next
	!assign id
        mTracer(kk+jj+1)%p%tracer_id=kk+jj+1	  
	!sever link between last tracer in current cell and this newly allocated head tracer
        nullify(mTracer(kk+jj)%p%next)
	!print*,'Tracer ',k+j,' link status:  ',associated(Tracer(k+j)%p%next)		
        kk=kk+jj+1			  
      end if
     end do
   else if(N_cell(ixx,iyy)==1)then
     !allocate memory for the head tracer of the next cell
     allocate(mTracer(kk)%p%next)
     mTracer(kk+1)%p=>mTracer(kk)%p%next
     !assign id
     mTracer(kk+1)%p%tracer_id=kk+1	  
     !sever link between last tracer in current cell and this newly allocated head tracer
     nullify(mTracer(kk)%p%next)
     kk=kk+1
   end if	  
  end if
 end do
end do

N_cell_temp= N_cell


!store initial tracer density distribution in file
!call tracerDensity()

print*,'Completed mass tracer initialization...' 

end subroutine massTracerInit

!-------------------------------------------------------------
subroutine tracerDensity(xmin,ymin,dx,dy)
!-------------------------------------------------------------
!store tracer distribution to file

!Input Variables
real*8::xmin,ymin,dx,dy

!Local variables
integer ::ix,iy
real*8::x,y
	
do ix=1,nx
 do iy=1,ny
  x=xmin+(ix-0.5)*dx
  y=ymin+(iy-0.5)*dy
  !write(13,*) x,y,N_cell(ix,iy)
 end do
end do

end subroutine tracerDensity


!-------------------------------------------------------------------------
subroutine massTracerAdvect(ix,iy,fL,fR,rho,dt,dr,dir) 
!-------------------------------------------------------------------------
!Input Variables
integer::ix,iy,dir
real*8::fL,fR,dt,dr,rho

!Local variables
integer::ixL=0,ixR=0,iyL=0,iyR=0 !destination cells
integer::N_out_L,N_out_R
real*8::m_cell,p_out,p_out1,rand_num,N_out_temp,N_p,N_p1,N_out_temp1

!print*,'Cell: ',ix,iy,', dir=',dir 
!print*,'No. of tracers at beginning of time step:  ',N_cell(ix,iy)


m_cell=rho*dr
p_out=rand(0)

if(N_cell(ix,iy)==0)then
N_out_L=0
N_out_R=0
else
 !compute number of tracers leaving (ix,iy)th cell 
 if(fL<0._8 .and. fR<0._8)then
    N_out_temp=(abs(fL)*dt/m_cell)*real(N_cell(ix,iy))
    N_p=N_out_temp-real(floor(N_out_temp))
    p_out=rand(0)
    if(p_out<N_p)then
      N_out_L=floor(N_out_temp)+1
    else
      N_out_L=floor(N_out_temp)
    end if
    N_out_R=0
    !print*,'m_cell,dt,N_out_L=',m_cell,dt,N_out_L   
! else
!   N_out_L =0
! end if
 else if(fR>0._8 .and. fL> 0._8)then
    N_out_temp=(abs(fR)*dt/m_cell)*real(N_cell(ix,iy))
    N_p=N_out_temp-real(floor(N_out_temp))
    p_out=rand(0)
    if(p_out<N_p)then
      N_out_R=floor(N_out_temp)+1
    else
      N_out_R=floor(N_out_temp)
    end if
    N_out_L=0
! else 
!  N_out_R =0
! end if
 !***Special Case*** Diverging flow in the cell 
 !(i.e. non-zero flux out of both cell faces), 
 else if(fL<0._8 .and. fR>0._8)then
    N_out_temp=(abs(fL)*dt/m_cell)*real(N_cell(ix,iy))
    N_p=N_out_temp-real(floor(N_out_temp))
    N_out_temp1=(abs(fR)*dt/m_cell)*real(N_cell(ix,iy))
    N_p1=N_out_temp-real(floor(N_out_temp1))
    p_out=rand(0)
    if(p_out<N_p .and. abs(fL)>abs(fR))then
      N_out_L=floor(N_out_temp)+1
    else
      N_out_L=floor(N_out_temp)
    end if
    if(p_out<N_p1 .and. abs(fR)>abs(fL))then
      N_out_R=floor(N_out_temp1)+1
    else
      N_out_R=floor(N_out_temp1)
    end if
 else
   N_out_L=0
   N_out_R=0
 end if

end if

		
if(trdebug==1 )then
  print*,'Cell: ',ix,iy,', N_out_L,N_out_R=',N_out_L,N_out_R 
  print*,'No. of tracers at beginning of time step:  ',N_cell(ix,iy)
  print*,'flux_left: ',fL,' , flux_right: ',fR
end if

!set destination cell
if(dir==1)then
 if(N_out_R>0) then
  if(trboundaryType==1 .or. trboundaryType==3)then
   if(ix==nx)then
     ixR=1
   else 
     ixR=ix+1
   end if
  else if(trboundaryType==2)then
    ixR=ix+1
  end if
 end if
 if(N_out_L>0) then
  if(trboundaryType==1 .or. trboundaryType==3)then
   if(ix==1)then
     ixL=nx
   else
     ixL=ix-1
   end if
  else if(trboundaryType==2)then
    ixL=ix-1
  end if
 end if
end if

if(dir==2)then
 if(N_out_R>0) then
  if(trboundaryType==1)then
   if(iy==ny)then
     iyR=1
   else 
     iyR=iy+1
   end if
  else if(trboundaryType==2 .or. trboundaryType==3)then
    iyR=iy+1
  end if
 else if(N_out_L>0) then
  if(trboundaryType==1)then
   if(iy==1)then
     iyL=ny
   else
     iyL=iy-1
   end if
  else if(trboundaryType==2 .or. trboundaryType==3)then
    iyL=iy-1
  end if
 end if
end if

!Update number of tracers remaining in relevant cells
N_cell_temp(ix,iy)=N_cell_temp(ix,iy)-N_out_L-N_out_R


!if(N_out_R>0 .and. N_out_L>0)then
!  print*,'Checkpoint 1'
!end if
	
!transfer tracers out of (ix,iy)th cell
if(dir==1)then 
 if(N_out_R>0) then !for positive flux, move tracers to right neighboring cell
   !print*,'Current Cell, Destination Cell:',ix,iy,',',ixR,iy,', Dir=',dir
   call cellTransfer(ix,iy,ixR,iy,N_out_R)  
 end if	
 N_cell(ix,iy)=N_cell(ix,iy)-N_out_R
 if(N_out_L>0) then !for negative flux, move tracers to left neighboring cell  
   !print*,'Current Cell, Destination Cell:',ix,iy,',',ixL,iy,', Dir=',dir
    call cellTransfer(ix,iy,ixL,iy,N_out_L)   
 end if	
else if(dir==2)then
 if(N_out_R>0) then !for positive flux, move tracers to top neighboring cell
   !print*,'Current Cell, Destination Cell:',ix,iy,',',ix,iyR,', Dir=',dir
   call cellTransfer(ix,iy,ix,iyR,N_out_R)  
 end if	
 N_cell(ix,iy)=N_cell(ix,iy)-N_out_R
 if(N_out_L>0) then !for negative flux, move tracers to bottom neighboring cell
   !print*,'Current Cell, Destination Cell:',ix,iy,',',ix,iyL,', Dir=',dir
   call cellTransfer(ix,iy,ix,iyL,N_out_L)   
 end if
end if


!if(N_out_R>0 .and. N_out_L>0)then
  !print*,'Checkpoint 2'
!end if

if((ix==nx .and. dir==1) .or. (iy==ny .and. dir==2))then
  N_cell=N_cell_temp
  !call tracerDensity()
end if

end subroutine massTracerAdvect

!-------------------------------------------------------------
subroutine cellTransfer(ix0,iy0,ix1,iy1,N_out)
!-------------------------------------------------------------
!Input Variables
integer ::ix0,iy0,ix1,iy1,N_out

!Local variables
integer::kk,ll,tr_counter,tr_counter2

!N_out tracers from the cell_(ix0,iy0) link are going to be moved into cell_(ix1,iy1).
!Pick a tracer from the link at random. Link it with the last tracer of the next cell.
!Then sever it's link with the tracer in the previous cell.
!Close the gap between tracers in previous cell that were adjacent to the tracer that's being removed.
!This is equivalent to moving the tracer from one cell to the next.
!Repeat this N_out times until the desired no. of tracers have been moved to the next cell. 
	
if(ix1>nx+1 .or. ix1<0 .or. iy1>ny+1 .or. iy1<0) then
 print*,'Cell: ',ix1,iy1
 print*,'Invalid destination cell index!'
 stop
end if
     		 
 do kk=1,N_out
   if(trdebug==1)then
     print*,'Number of tracers left to move out:',N_out-kk+1     
     print*,'Current Cell, Destination Cell:',ix0,iy0,',',ix1,iy1 
   end if
	  
   tr_counter=0  
   !traverse through list until last tracer is reached
   tracer_temp=>cell_tracer_head(ix0,iy0)%p
   if(trdebug==1)then	
     print*,'tracer_temp currently pointing to tracer#',tracer_temp%tracer_id
   end if
   !print*,'tracer_temp inter=',N_cell(ii)-kk
   do while(tr_counter<N_cell(ix0,iy0)-kk)! .and. associated(tracer_temp%next))
     tracer_temp=>tracer_temp%next
     tr_counter=tr_counter+1
     if(trdebug==1)then
       print*,'N_cell=',N_cell(ix0,iy0)
       print*,'tr_counter=',tr_counter 	 
       print*,'tracer_temp currently pointing to tracer#',tracer_temp%tracer_id
     end if		  
   end do !loop ends when tracer_temp is pointing to the last tracer in link
	  !tr_counter holds position of last tracer in the link
	  !(i.e. # of tracers in cell)

   !generate a random number to decide which tracer in the link will be moved		
   tr_counter2= tr_counter*rand(0)+1

   if(trdebug==1)then
     print*,'Tracer#',tr_counter2,'in list is going to be moved.'
   end if		
   !need to be careful if moving the head tracer in the list

   !moving a tracer that is not the head or last tracer
   if(tr_counter2>1 .and. tr_counter2<=tr_counter+1) then 
     !traverse through list until this tracer is reached
     tracer_temp=>cell_tracer_head(ix0,iy0)%p
     do ll=1,tr_counter2-1
       tracer_temp=>tracer_temp%next		  
     end do	
     tracer_temp3=>tracer_temp%next     
     !now traverse to the last tracer of the destination cell
     tracer_temp2=>cell_tracer_head(ix1,iy1)%p

     !Update number of tracers in next cell
     N_cell_temp(ix1,iy1)=N_cell_temp(ix1,iy1)+1
     !record new location of the tracer that's been moved
     tracer_temp%cell_occupied(1)=ix1
     tracer_temp%cell_occupied(2)=iy1

     if(trdebug==1) then
       print*,'(Case 1)Moving tracer#', tracer_temp%tracer_id, 'from cell ', &
       ix0,iy0,' to cell ',ix1,iy1,'tr_counter= ',tr_counter
     end if		
  
     !if next cell is empty then link its head tracer with the tracer being moved
     if(.not. associated(cell_tracer_head(ix1,iy1)%p))then
       if(trdebug==1)then
         print*,'Cell#',ix1,iy1,'is empty...'
       end if
       cell_tracer_head(ix1,iy1)%p=>tracer_temp
     !if next cell not empty, then link last tracer in that cell with the tracer being moved
     else
       do while(associated(tracer_temp2%next))
         tracer_temp2=>tracer_temp2%next
       end do !loop ends when tracer_temp2 is pointing to the last tracer in link

       !link last tracer of next cell with last tracer of previous cell
       allocate(tracer_temp2%next)
       tracer_temp2%next=>tracer_temp 
     end if
		
     !sever link with previous cell
     tracer_temp=>cell_tracer_head(ix0,iy0)%p		 
     do ll=1,tr_counter2-1      
       tracer_temp=>tracer_temp%next  
     end do
     nullify(tracer_temp%next)
		
     !Close the gap between tracers in previous cell that were adjacent to the tracer that's being removed
     tracer_temp=>cell_tracer_head(ix0,iy0)%p	
     do ll=1,tr_counter2-2       
       tracer_temp=>tracer_temp%next  
     end do
     tracer_temp%next=>tracer_temp3
     !print*,'Linking tracer#',tracer_temp%tracer_id,' with tracer#',tracer_temp3%tracer_id		  
   end if
	
   !moving the head tracer
   if(tr_counter2==1) then 
     !print*,'***Case 2***'  
     !point to head tracer which is going to be moved
     tracer_temp=>cell_tracer_head(ix0,iy0)%p
     !promote the next tracer in link to new head tracer
     cell_tracer_head(ix0,iy0)%p=>tracer_temp%next
     !print*,'Promoting tracer#',cell_tracer_head(i)%p%tracer_id, ' to new head tracer for cell#',i	
			
     !now traverse to the last tracer of the next cell
     tracer_temp2=>cell_tracer_head(ix1,iy1)%p
     !Update number of tracers in next cell
     N_cell_temp(ix1,iy1)=N_cell_temp(ix1,iy1)+1
     !record new location of the tracer that's been moved
     tracer_temp%cell_occupied(1)=ix1
     tracer_temp%cell_occupied(2)=iy1
 
     if(trdebug==1)then		 
       print*,'(Case 2)Moving tracer#', tracer_temp%tracer_id, 'from cell ', &
       ix0,iy0,' to cell ',ix1,iy1,'tr_counter= ',tr_counter 
     end if
        
     !if next cell is empty then link head tracer head tracer with the tracer being moved
     if(.not. associated(cell_tracer_head(ix1,iy1)%p))then
     !if(N_cell_temp(ix1,iy1)-1==0)then
       if(trdebug==1)then
         print*,'Cell#',ix1,iy1,'is empty...'
       end if
       cell_tracer_head(ix1,iy1)%p=>tracer_temp
       !if next cell not empty, then link last tracer in that cell with the tracer being moved
     else 
       do while(associated(tracer_temp2%next))
         tracer_temp2=>tracer_temp2%next
       end do !loop ends when tracer_temp2 is pointing to the last tracer in link
				 
       !link last tracer of next cell with last tracer of previous cell
       allocate(tracer_temp2%next)
     	 tracer_temp2%next=>tracer_temp 
     end if
     nullify(tracer_temp%next)		
   end if
 end do	

end subroutine cellTransfer

subroutine massTracerOutput(fileunit)
!Input variables
integer::fileunit
!Local variables
integer::jj,cell(3)
real*8::vel(3),dens,div

do jj=1,N
 cell(1)=mTracer(jj)%p%cell_occupied(1)
 cell(2)=mTracer(jj)%p%cell_occupied(2)
 vel(1)=u1(cell(1),cell(2),2)/u1(cell(1),cell(2),1)
 vel(2)=u1(cell(1),cell(2),3)/u1(cell(1),cell(2),1)
 dens=u1(cell(1),cell(2),1)
 div=divU(cell(1),cell(2))
 write(fileunit,*) mTracer(jj)%p%tracer_id,cell(1), &
                   cell(2),vel(1),vel(2),dens,div

end do

end subroutine massTracerOutput

end module mass_tracer_mod_2D
