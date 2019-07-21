!Chang-Cooper(1970) algorithm for Isotropic Cosmic Ray-Fokker Planck Equation
!This code evolves CR momentum space evolution on individual tracer particles 
!i.e. on a lagrangian spatial frame 

!Note1: Momentum (p) is in units of m_p c (where m_p is the CR particle mass 
!and c is the speed of light)


module FokkerPlanck_mod
use constants_mod
use readTracer_module

implicit none

integer,parameter::subcycleOption=0 !0: off, 1:on
real*8::dt2,dt3,dtf
real*8,parameter::cour2=0.5
integer::nt2

contains

subroutine CRinit()

!Local variables
integer::i,j

open(unit=20,file='f(p)_1.txt')

!----------------------------------------------
!Initialization
!----------------------------------------------
!Momentum space bounds
pmin=1.D-2
pmax=1.D10

print*,'np,pmin,pmax=',np,pmin,pmax

nparticle=0.

!Generate momentum mesh points
do i=-2,np+3
  if(momentumMeshType==1)then
   px(i)=pmin+i*(pmax-pmin)/np
  else if(momentumMeshType==2)then
   px(i)=pmin*(pmax/pmin)**(real(i)/real(np-1))
  end if
end do


!Momentum mesh spacing
do i=-1,np+2
 dp(i)=0.5*(px(i+1)-px(i-1))
 dp2(i)=px(i)-px(i-1)
 !print*,'i,p,dp=',i,px(i),dp(i)
end do

!Check diffusion length and spatial resolution requirement
!xd=Kxx(0,0)
!do i=0,np-1
! do j=0,nx-1
!  xd=min(xd,Kxx(j,i))
! end do
!end do
!print*,'xd,dx',xd,dx!,dx/xd

do i=1,N
 do j=0,np-1
  !Initial momentum distribution
  !Power law
  f(i,j)=px(j)**(-4.5)
  !Gaussian/Thermal
  !f(i,j)=1000.*exp(-( (pmin+j*dp(j)-0.5*(pmax-pmin))/(20.*dp(j)) )**2. ) 
 end do
 f(i,-1)=f(i,0)  
 f(i,np)=f(i,np-1)
end do

do j=0,np-1
  write(20,*) pmin+j*dp(j),f(1,j) 
end do



print*,'Done initializing CR distribution.'

end subroutine

subroutine CRevolve()
!Evolve Distribution Function
integer::i,j,k
real*8::ftemp(-2:np+1)

!-------------------------------------------------------
!Mometum Operator Sweeps (with sub-cycled time stepping)   
!-------------------------------------------------------
 if(subcycleOption==1)then
 !Momentum advection time step
 dt2=dt*2. 
 
 do j=1,N
  if(abs(file_dat%divU(j))>1.d-20)then
   do k=0,np-1
    dt2=min(dt2,cour2*3.*dp(k)/(px(k)*abs(file_dat%divU(j))))
   end do
  else
   dt2=min(dt2,dt)
  end if
 end do

 if(dt2<dt)then
  nt2=floor(dt/dt2)
  dtf=dt-nt2*dt2
  print*,'Momentum advection time step, dt2=',dt2
  print*,'# of momentum sweep sub cycles=',nt2+1
 else
   nt2=1
 end if
    
 dt3=dt2

 else if(subcycleOption==0)then
  nt2=0
  dt2=dt
  dt3=dt
  dtf=0.
 end if


do i=1,nt2+1

  do j=1,N
   ftemp=0.
   do k=-1,np
    ftemp(k)=f(j,k)
   end do
   

   call computeMomentumCoefficients(j)

   call tridiagSolveMomentum(np,ftemp,j)


   ftemp(-1)=ftemp(0)  
   ftemp(np)=ftemp(np-1)

   do k=-1,np
    f(j,k)=ftemp(k)  
   end do

  end do
  if(i==nt2) dt2=dtf
end do

!----------------------------------------------
!Spatial Operator Sweeps: Spatial diffusion only   
!No spatial advection, assuming no CR particle
!drift motion relative to the lagrangian frame  
!----------------------------------------------

!spatial diffusion, coming soon...

  
!Save momentum distribution to file

!Tracer 1
do j=0,np-1
  write(20,*) pmin+j*dp(j),f(1,j) 
end do



end subroutine CRevolve

subroutine computeMomentumCoefficients(xk)

!Input variables
integer::xk

integer::jj
real*8::w(-2:np+2),Wplus(-1:np+1),Wminus(-1:np+1),Wbig(-1:np+1)
real*8::delta(-2:np+2),Cp(-1:np+1)
real*8::B2,C2,B1,C1,pi2,pi2L,w1,w2

!-------------------------------------------------------------
!Weights: w_j, delta_j, W^+_j+1/2, W^-_j+1/2 
!-------------------------------------------------------------
do jj=-1,np+1
 C2=0.5*(C(jj)+C(jj+1))
 B2=0.5*(B(jj,xk)+B(jj+1,xk))


 if(abs(C2)<1.d-20)then
  if(B2<0._8)then
   delta(jj)=1.
  else
   delta(jj)=0._8
  end if
 else
  w2=dp2(jj)*B2/C2
  if(abs(w2)<1.d-3)then
    delta(jj)=0.5-(1./12.)*w2   
  else if(abs(w2)>=700.)then
    delta(jj)=(1./w2) !to prevent floating-point overflow
  else  
    !print*,'1/w2,exp(w2)-1.=',(1./w2),(exp(w2)-1.)
    delta(jj)=(1./w2)-(1./(exp(w2)-1.))
  end if

 end if

end do

do jj=-1,np+1
  B2=0.5*(B(jj,xk)+B(jj+1,xk))
  Cp(jj)=0.5*(C(jj)+C(jj+1))/dp2(jj)
  Wminus(jj)=B2*delta(jj)
  Wplus(jj)=B2*(1.-delta(jj))
end do

!-------------------------------------
!No-flux/Continuous Boundary condition 
!-------------------------------------

!for no-flux boundary on the left
Wminus(-1)=0.
Wplus(-1)=0.
Cp(-1)=0.

!for no-flux boundary on the right
Wplus(np-1)=0.
Wminus(np-1)=0.
Cp(np-1)=0.

!------------------------------------------------------------
!Co-efficients for tridiagonal system: ~A_j,~B_j,~C_j
!------------------------------------------------------------
do jj=0,np-1

 At(jj)=dt2/(dp(jj)*A(jj))
 At(jj)=At(jj)*(Wminus(jj-1)-Cp(jj-1))
 if(methodType==2)then
   At(jj)=At(jj)*0.5 
  end if

 Ct(jj)=dt2/(dp(jj)*A(jj))
 Ct(jj)=-Ct(jj)*(Wplus(jj)+Cp(jj))
 if(methodType==2)then
   Bt(jj)=Bt(jj)*0.5 
 end if  

 Bt(jj)= Wminus(jj)-Wplus(jj-1)-Cp(jj)-Cp(jj-1)
 Bt(jj)=Bt(jj)*dt2/(dp(jj)*A(jj))
 Bt(jj)=-Bt(jj)+1.-dt2*Hp(xk)
 if(methodType==2)then 
   Bt(jj)=0.5*(Bt(jj)+1.)
 end if   

end do
 
end subroutine computeMomentumCoefficients

subroutine tridiagSolveMomentum(n,ftemp,xk) 
!Gaussian elimination algorithm fomr Numerical Recipes book

!Input variables
integer::n,xk
real*8::ftemp(-1:n)

real*8::bet,gam(0:n-1),r(0:n-1)

integer::kk


do kk=0,n-1
 r(kk)=ftemp(kk)+dt2*S(kk,xk)
 if(methodType==2)then
   r(kk)=r(kk)+ftemp(kk)-At(kk)*ftemp(kk-1)-Bt(kk)*ftemp(kk)-Ct(kk)*ftemp(kk+1)
 end if
end do


bet=Bt(0)
ftemp(0)=r(0)/bet

!Upper triangular decomposition and forward substitution
do kk=1,n-1
  gam(kk)=Ct(kk-1)/bet
  bet=Bt(kk)-At(kk)*gam(kk)
  ftemp(kk)=(r(kk)-At(kk)*ftemp(kk-1))/bet
end do

!Back substitution
do kk=n-2,0,-1
 ftemp(kk)=ftemp(kk)-gam(kk+1)*ftemp(kk+1)
end do

end subroutine tridiagSolveMomentum



!Momentum Phase factor
function A(pk) result(ax)
integer,intent(in)::pk
real*8::ax
  
 ax=px(pk)*px(pk)

end function A

!Momentum Advection co-efficient
function B(pk,xk) result(bx)
integer,intent(in)::pk,xk
real*8::bx

bx=(px(pk)**(3.))/3.
bx=bx*file_dat%divU(xk)

end function B

!Momentum Diffusion co-efficient
function C(pk) result(cx)

integer,intent(in)::pk
real*8::cx,Dpp

 Dpp=1.0
 cx=Dpp*px(pk)*px(pk)

end function C

!Spatial Diffusion co-efficient
function Kxx(xk,pk) result(kx)

integer,intent(in)::xk,pk
real*8::kx

kx=0.05

end function Kxx


!Momentum Source term co-efficient
function Hp(xk) result(hx)
 
integer,intent(in)::xk
real*8::hx

hx=-divU(xk)

end function Hp


!Spatial Source term co-efficient
function Hx(xk) result(hxx)
 
integer,intent(in)::xk
real*8::hxx

hxx=0.!-dudx(xk)/3.

end function Hx


!Momentum Space Injection term
function S(pk,xk) result(sx)
 
integer,intent(in)::pk,xk
real*8::sx

real*8,parameter::alpha=2.
real*8::eps=0.005
real*8::mp=0.1

real*8::cs2,rho1,rho2,pinj,us,pres,wx
real*8::ut1,ut2
integer::pin

!------------------------------
!Flux fraction injection model
!------------------------------
!Shock Speed
!ut1=-10.*sqrt(gam)    !u1(shockCell(1),2)/u1(shockCell(1),1)
!ut2=0.               !u1(shockCell(1)-2,2)/u1(shockCell(1)-2,1)
!us=(1./3.)*(4.*ut2-ut1)     !(1./3.)*2.*sqrt(gam)
!pre-shock density
!rho1=1.!u1(shockCell(1),1)
!Post-shock sound speed
!pres=(gam-1.)*(u1(shockCell(1)-1,3) &
!-0.5*u1(shockCell(1)-1,2)*u1(shockCell(1)-1,2)/u1(shockCell(1)-1,1))
!rho2=u1(shockCell(1)-1,1)
!cs2=sqrt(gam*pres/rho2)

!Compute injection momentum

!pinj=alpha*cs2
!pin=np*log(pinj/pmin)/log(pmax/pmin)

!print*,'pinj,pin=',pinj,pin


!sx=0.
!injection at shock location
!if(xk==shockCell(1)-1)then
!  sx=(px(pk)-100.)/(10.)
!  sx=exp(-sx**2.)

!if(pk==pin .and. pin>0 .and. pin <np-1 )then
 
 !Spatial weight function
! wx=(dx*xk-dx*(shockCell(1)-1))/(4.*dx)
! wx=exp(-wx**2.)/sqrt(3.141592*(4.*dx)**2.)
 !if(xk==shockCell(1)-1)then
 !  wx=1.
 !else
 !  wx=0.
 !end if
 !Injection source term
! sx=wx*eps*rho1*us/(4.*3.141592*px(pk)*px(pk))/mp 
!else
!  sx=0._8
!end if

sx=0._8

end function S

!Spatial Injection term
function L(xk,t) result(lx)
 
integer,intent(in)::xk
real*8,intent(in)::t
real*8::lx

lx=0.

end function L

end module FokkerPlanck_mod
