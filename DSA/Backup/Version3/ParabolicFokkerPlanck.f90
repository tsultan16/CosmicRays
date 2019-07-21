!Chang-Cooper(1970) algorithm for Isotropic Cosmic Ray-Fokker Planck Equation

module FokkerPlanck_mod
use constants_mod
implicit none

integer,parameter::subcycleOption=0 !0: off, 1:on
real*8::fexact(0:np-1),dt2,dt3,dtf
real*8,parameter::cour2=0.5
integer::nt2

contains

subroutine CRinit()

integer::i,j
real*8::pres

open(unit=10,file='fp.txt')
open(unit=13,file='fx.txt')
!----------------------------------------------
!Initialization
!----------------------------------------------
!Momentum space bounds
pmin=1.D0
pmax=1.D20

print*,'np,pmin,pmax=',np,pmin,pmax

nparticle=0.

!Generate momentum mesh points
do i=-2,np+3
  if(momentumMeshType==1)then
   px(i)=pmin+i*(pmax-pmin)/np
  else if(momentumMeshType==2)then
   px(i)=pmin*(pmax/pmin)**(real(i)/real(np))
  end if
  !print*,'i,px(i)=',i,px(i)
end do


!Momentum mesh spacing
do i=-1,np+2
 dp(i)=0.5*(px(i+1)-px(i-1))
 dp2(i)=px(i+1)-px(i)
 !print*,'i,p,dp=',i,px(i),dp(i)
end do


do i=0,np-1
 do j=0,nx-1
  !Initial profile
  !if(px(i)>=0.1)then
  ! f(i,j)=px(i)**(-6.)!*exp(-((xmin+j*dx)/(10.*dx))**2. )
  !end if
  f(i,j)=0._8
  !f(i,j)=1000.*exp(-( (xmin+j*dx-0.5*(xmax-xmin))/(5.*dx) )**2. )!
  !f(i,j)=1000.*exp(-( (pmin+i*dp(i)-0.5*(pmax-pmin))/(20.*dp(i)) )**2. ) 
 end do
 f(-1,j)=f(0,j)  
 f(np,j)=f(np-1,j)
end do



!Normalize Distribution function such that CR pressure is equal to downstream fluid pressure
do j=0,nx-1
 do i=0,np-1
  Pc(j)=Pc(j)+(4.*cl/3.)*(px(i)**3.)*f(i,j)*dp(i)
 end do
end do

pres=(gam-1.)*(u1(nx,3)-0.5*u1(nx,2)*u1(nx,2)/u1(nx,1))

do i=0,np-1
 do j=0,nx-1
  !f(i,j)=f(i,j)*(pres/Pc(j))
  fexact(i)=4.*f(i,j)
 end do
end do


fp=0.
do i=0,np-1
 do j=0,nx-1
  fp(i)=fp(i)+dx*f(i,j)
 end do
 write(10,*) px(i),fp(i)+1.D-50,fexact(i)
 !print*,'p,f(p)',px(i),fp(i)
end do

fx=0.
do j=0,nx-1
  do i=0,np-1
  fx(j)=fx(j)+dp(i)*f(i,j)
 end do
 write(13,*) xmin+j*dx,fx(j)!f(1,j)!fx(j)!+1.D-50
 !print*,'x,f(x)',xmin+j*dx,fx(j)
end do


!Compute CR Pressure
Pc=0.
do j=0,nx-1
  do i=0,np-1
    Pc(j)=Pc(j)+(4.*cl/3.)*(px(i)**3.)*f(i,j)*dp(i)
  end do
end do
Pc(-1)=Pc(0)
Pc(nx)=Pc(nx-1)


!do j=0,nx-1
!   pres=(gam-1.)*(u1(j,3)-0.5*u1(j,2)*u1(j,2)/u1(j,1))
!  print*,'j,Pc,Pg=',j,Pc(j),pres
!end do

print*,'Done initializing CR distribution.'

end subroutine

subroutine CRevolve()
!Evolve Distribution Function
integer::i,j,k

!-------------------------------------------------------
!Mometum Operator Sweeps (with sub-cycled time stepping)   
!-------------------------------------------------------
 if(subcycleOption==1)then
 !Momentum advection time step
 dt2=dt*2. 
 
 do j=0,nx-1
  if(abs(dudx(j))>1.d-20)then
   do k=0,np-1
    dt2=min(dt2,cour2*3.*dp(k)/(px(k)*abs(dudx(j))))
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
 end if

 do i=1,nt2+1

  do j=0,nx-1
   ftemp1=0.
   do k=-1,np
    ftemp1(k)=f(k,j)
   end do
   
   call computeMomentumCoefficients(j)
   call tridiagSolveMomentum(np,ftemp1,j)
   
   ftemp1(-1)=ftemp1(0)  
   ftemp1(np)=ftemp1(np-1)

   do k=-1,np
    f(k,j)=ftemp1(k)  
    !print*,'j,k,ftemp1=',j,k,ftemp1(k) 
   end do

  end do

  if(i==nt2) dt2=dtf
 end do

!print*,'CHECKPOINT1'
!-------------------------------------
!Spatial Operator Sweeps   
!-------------------------------------
do i=1,nt2+1
 
 do j=0,np-1
  ftemp2=0.
  do k=-1,nx
   ftemp2(k)=f(j,k)
  end do  

  !print*,'CHECKPOINT1'
  call computeSpatialCoefficients(j)
  !print*,'CHECKPOINT2'
  call tridiagSolveSpace(nx,ftemp2)

  ftemp2(-1)=ftemp2(0)  
  ftemp2(nx)=ftemp2(nx-1)

  do k=-1,nx
  f(j,k)=ftemp2(k)   
  end do
 end do

 if(i==nt2) dt3=dtf
end do
!print*,'CHECKPOINT2'

!-----------------------------------------------
!Update CR pressure
!-----------------------------------------------
Pc=0.
do j=0,nx-1
  do k=0,np-1
    Pc(j)=Pc(j)+(4.*cl/3.)*(px(k)**3.)*f(k,j)*dp(k)
  end do
end do
Pc(-1)=Pc(0)
Pc(nx)=Pc(nx-1)


fp=0.
do k=0,np-1
 do j=0,nx-1
  fp(k)=fp(k)+dx*f(k,j)
 end do
end do

fx=0.
do j=0,nx-1
 do k=-1,np
  fx(j)=fx(j)+dp(k)*f(k,j)
 end do
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
if(methodType .ne. 4)then
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
  else
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

!---------------------------------
!No-flux Boundary condition 
!---------------------------------
Wplus(0)=0.
Wminus(0)=0.
Cp(0)=0.
Wplus(-1)=0.
Wminus(-1)=0.
Cp(-1)=0.
Wplus(np-1)=0.
Wminus(np-1)=0.
Cp(np-1)=0.

end if

!------------------------------------------------------------
!Co-efficients for tridiagonal system: ~A_j,~B_j,~C_j
!------------------------------------------------------------
do jj=0,np-1
  if(methodType .ne. 4)then
   At(jj)=dt2/(dp(jj)*A(jj))
   At(jj)=At(jj)*(Wminus(jj-1)-Cp(jj-1))
  end if

  if(methodType==4)then
  At(jj)=dt2/(dp(jj)*dp2(jj-1)*A(jj))
  At(jj)=At(jj)*(C1-0.5*B1*dp2(jj-1))
  At(jj)=-At(jj)
  end if

  if(methodType .ne. 4)then
   Ct(jj)=dt2/(dp(jj)*A(jj))
   Ct(jj)=-Ct(jj)*(Wplus(jj)+Cp(jj))
  end if
  if(methodType==4)then
   Ct(jj)=dt2/(dp(jj)*dp2(jj)*A(jj))
   Ct(jj)=Ct(jj)*(C2+0.5*B2*dp2(jj))
   Ct(jj)=-Ct(jj)
  end if

  if(methodType .ne. 4)then
   Bt(jj)= Wminus(jj)-Wplus(jj-1)-Cp(jj)-Cp(jj-1)
   Bt(jj)=Bt(jj)*dt2/(dp(jj)*A(jj))
   Bt(jj)=-Bt(jj)+1.-dt2*Hp(xk)
  end if
  if(methodType==4)then
   Bt(jj)=(C1+0.5*B1*dp2(jj-1))/dp2(jj-1) &
        +(C2+0.5*B2*dp2(jj))/dp2(jj)
   Bt(jj)=Bt(jj)*dt2/(dp(jj)*A(jj))
   Bt(jj)=Bt(jj)+1.+dt2*Hp(xk)
  end if 


  !print*,'jj,At,Bt,Ct=',jj,At(jj),Bt(jj),Ct(jj)
end do
 
end subroutine computeMomentumCoefficients

subroutine tridiagSolveMomentum(n,ftemp,xk) !from Numerical Recipes book

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


subroutine computeSpatialCoefficients(pk)

!Input variables
integer::pk,xk
integer::jj
real*8::Wplus(-1:nx),Wminus(-1:nx)
real*8::B2,C2,B1,C1,pi2,pi2L,w2
real*8::delta(-1:nx),Cp(-1:nx)

!-------------------------------------------------------------
!Weights: w_j, delta_j, W^+_j+1/2, W^-_j+1/2 
!-------------------------------------------------------------
do jj=-1,nx
 B2=0.5*(E(jj)+E(jj+1))
 C2=Kxx(jj,pk)

 !print*,'jj,B2,C2=',jj,B2,C2

 if(abs(C2)<1.d-20)then
  if(B2<0._8)then
   delta(jj)=1.
  else
   delta(jj)=0._8
  end if
 else
  w2=dx*B2/C2
  if(abs(w2)<1.d-3)then
    delta(jj)=0.5-(1./12.)*w2   
  else
    delta(jj)=(1./w2)-(1./(exp(w2)-1.))
  end if
 end if

end do

do jj=-1,nx
  B2=0.5*(E(jj)+E(jj+1))
  Cp(jj)=Kxx(jj,pk)/dx
  Wminus(jj)=B2*delta(jj)
  Wplus(jj)=B2*(1.-delta(jj))
end do

!---------------------------------
!No-flux Boundary condition 
!---------------------------------

Wplus(0)=0.
Wminus(0)=0.
Cp(0)=0.
Wplus(-1)=0.
Wminus(-1)=0.
Cp(-1)=0.
!Wplus(nx-1)=0.
!Wminus(nx-1)=0.
!Cp(nx-1)=0.



!------------------------------------------------------------
!Co-efficients for tridiagonal system: ~A_j,~B_j,~C_j
!------------------------------------------------------------
do jj=0,nx-1
  Atx(jj)=dt3/dx
  Atx(jj)=Atx(jj)*(Wminus(jj-1)-Cp(jj-1))
  if(methodType==2)then
   Atx(jj)=Atx(jj)*0.5 
  end if
  
  Ctx(jj)=dt3/dx
  Ctx(jj)=-Ctx(jj)*(Wplus(jj)+Cp(jj))
  if(methodType==2)then
   Ctx(jj)=Ctx(jj)*0.5
  end if

   Btx(jj)= Wminus(jj)-Wplus(jj-1)-Cp(jj)-Cp(jj-1)
   Btx(jj)=Btx(jj)*dt3/dx
   Btx(jj)=-Btx(jj)+1.-dt3*Hx(jj)

  if(methodType==2)then 
   Btx(jj)=0.5*(Btx(jj)+1.)
  end if

  !print*,'jj,Atx,Btx,Ctx=',jj,Atx(jj),Btx(jj),Ctx(jj)
end do


end subroutine computeSpatialCoefficients



subroutine tridiagSolveSpace(n,ftemp) !from Numerical Recipes book

!Input variables
integer::n
real*8::ftemp(-1:n)

real*8::bet,gam(0:n-1),r(0:n-1)

integer::kk


do kk=0,n-1
 !print*,'k,p,S(p,t)=',kk,px(kk),S(kk,t)
 r(kk)=ftemp(kk)+dt3*L(kk,t)
 if(methodType==2)then
   r(kk)=r(kk)+ftemp(kk)-Atx(kk)*ftemp(kk-1)-Btx(kk)*ftemp(kk)-Ctx(kk)*ftemp(kk+1)
 end if
end do


bet=Btx(0)
ftemp(0)=r(0)/bet

!Upper triangular decomposition and forward substitution
do kk=1,n-1
  gam(kk)=Ctx(kk-1)/bet
  bet=Btx(kk)-Atx(kk)*gam(kk)
  ftemp(kk)=(r(kk)-Atx(kk)*ftemp(kk-1))/bet
end do

!Back substitution
do kk=n-2,0,-1
 ftemp(kk)=ftemp(kk)-gam(kk+1)*ftemp(kk+1)
end do

end subroutine tridiagSolveSpace

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

bx=(px(pk)**(3.))*dudx(xk)/3.


end function B

!Momentum Diffusion co-efficient
function C(pk) result(cx)

integer,intent(in)::pk
real*8::cx,Dpp

 Dpp=0.
 cx=Dpp*px(pk)*px(pk)

end function C

!Spatial Advection co-efficient
function E(xk) result(ex)
integer,intent(in)::xk
real*8::ex

ex=-u1(xk,2)/u1(xk,1)

end function E

!Spatial Diffusion co-efficient
function Kxx(xk,pk) result(kx)

integer,intent(in)::xk,pk
real*8::kx

kx=0.5/u1(xk,1)
!kx=0.1*(px(pk)**0.51)


end function Kxx


!Momentum Source term co-efficient
function Hp(xk) result(hx)
 
integer,intent(in)::xk
real*8::hx

hx=0.

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

real*8,parameter::lambda=2.
real*8::eps=0.1
real*8::mp=0.1

real*8::pinj,us,pres,wx
integer::pin

!------------------------------
!Flux fraction injection model
!------------------------------
!Compute Shock Speed
us=gam**1.5

!Compute injection momentum
pres=55.78043!(gam-1.)*(u2(0,3)-0.5*u2(0,2)*u2(0,2)/u2(0,1))
pinj=lambda*sqrt(gam*pres/ 3.976468)
pin=np*log(pinj/pmin)/log(pmax/pmin)

!print*,'pinj,pin=',pinj,pin


!sx=0.
!injection at shock location
!if(xk==shockCell(1)-1)then
!  sx=(px(pk)-100.)/(10.)
!  sx=exp(-sx**2.)
if(pk==pin .and. pin>0 .and. pin <np-1 )then
 
 !Spatial weight function (delta function)
 
 !wx=(dx*xk-dx*(shockCell(1)-1))/(1.*dx)
 !wx=exp(-sx**2.)/sqrt(3.141592*sqrt(10.*dx))
 if(xk==shockCell(1)-1)then
   wx=1.
 else
   wx=0.
 end if
 !Injection source term
 sx=wx*eps*us/(4.*3.141592*px(pk)*px(pk))/mp 
else
  sx=0._8
end if



end function S

!Spatial Injection term
function L(xk,t) result(lx)
 
integer,intent(in)::xk
real*8,intent(in)::t
real*8::lx

lx=0.

end function L

function fluidS(xk) result(fx)
integer,intent(in)::xk
real*8::wx

real*8,parameter::lambda=2.
real*8::eps=0.1
real*8::mp=0.1

real*8::pinj,us,pres,fx
integer::pin

!------------------------------
!Flux fraction injection model
!------------------------------
!Compute Shock Speed
us=gam**1.5

!Compute injection momentum
pres=55.78043!(gam-1.)*(u2(0,3)-0.5*u2(0,2)*u2(0,2)/u2(0,1))
pinj=lambda*sqrt(gam*pres/ 3.976468)
pin=np*log(pinj/pmin)/log(pmax/pmin)

if(xk==shockCell(1)-1)then
   wx=1.
 else
   wx=0.
 end if

fx=0.5*eps*wx*(pin**2.)*us!/mp


end function fluidS

subroutine computePFlux(ftemp1,xk)
!Input Variables
integer::xk
real*8::ftemp1(-1:np)
!Local Variables
integer::jj,fupdate(-1:np),pflux(-1:np+2)
real*8::w(-2:np+2),Wplus(-1:np+1),Wminus(-1:np+1),Wbig(-1:np+1)
real*8::B2,C2,B1,C1,pi2,pi2L

!-------------------------------------------------------------
!Weights: w_j, W_j, W^+_j+1/2, W^-_j+1/2 , F^n_i+1/2
!-------------------------------------------------------------
do jj=-1,np+2
 B2=B(jj,xk)
 C2=C(jj)
 w(jj)=dp(jj)*B2/C2 

 if(abs(w(jj))>1.d4 .or. abs(w(jj))<1.d-12)then
   print*,'***ERROR***, w_j=',w(jj)
   print*,'dp,B2,C2=',dp(jj),B2,C2
   print*,'w_j is out of bounds.' 
   print*,'Need to re adjust size of momentum advection and diffusion co-efficients.' 
   STOP
 end if
 Wbig(jj)=0.5*w(jj)/sinh(0.5*w(jj))
end do

do jj=-1,np+1
  Wminus(jj)=0.5*(Wbig(jj)*exp(-0.5*w(jj))+Wbig(jj+1)*exp(-0.5*w(jj+1)) )
  Wplus(jj)=0.5*(Wbig(jj)*exp(0.5*w(jj))+Wbig(jj+1)*exp(0.5*w(jj+1)) )
end do

do jj=-1,np-1
 C2=C(jj)
 pflux(jj)=(C2/dp2(jj))*(WPlus(jj)*ftemp1(jj+1)-WMinus(jj)*ftemp1(jj))
end do

!No-flux Boundary Condition
pflux(-1)=0.
pflux(np-1)=0.

!Update distribution function
!fupdate=ftemp1
do jj=0,np-1
ftemp1(jj)=ftemp1(jj)*(1.+dt*Hp(jj)) &
           +(dt/A(jj))*(pflux(jj)-pflux(jj-1))/dp(jj) &
           + dt*S(jj,xk)
end do
!ftemp1=fupdate

end subroutine computePFlux

subroutine spaceAdvect(ftemp)
!Input Variables
real*8::ftemp(-1:nx)

real*8::ftemp1(-1:nx),a
integer::jj

ftemp1=ftemp

!Spatial advection: First Order Upwind
!fupdate=ftemp1
do jj=0,nx-1
!a=0.5*(u2(jj,2)/u2(jj,1)+u1(jj,2)/u1(jj,1))
a=u1(jj,2)/u1(jj,1)
! if(u1(jj,2)<0._8)then
!  ftemp(jj)=ftemp(jj)-(dt/dx)*(u1(jj,2)/u1(jj,1))*(ftemp1(jj+1)-ftemp1(jj))
! else
!  ftemp(jj)=ftemp(jj)-(dt/dx)*(u1(jj,2)/u1(jj,1))*(ftemp1(jj)-ftemp1(jj-1))
! end if
  !print*,'jj,a=',jj,a
  ftemp(jj)=ftemp(jj)-0.5*(dt/dx)*a*( &
           (1._8-sign(1._8,a))*(ftemp1(jj+1)-ftemp1(jj)) &
          +(1._8+sign(1._8,a))*(ftemp1(jj)-ftemp1(jj-1)) )
end do

!Zero flux out of left-boundary
!if(u1(0,2)<0._8)then
! ftemp(0)=ftemp(0)-(dt/dx)*(u1(0,2)/u1(0,1))*(ftemp1(0+1))
!else
! ftemp(0)=ftemp(0)-(dt/dx)*(u1(0,2)/u1(0,1))*(ftemp1(0))
!end if


end subroutine spaceAdvect

subroutine momentumAdvect(xk,ftemp)
!Input Variables
integer::xk
real*8::ftemp(-1:np)

!Local Variables
integer::jj
real*8::B2,ftemp1(-1:np),a,pflux(-1:np+2)


ftemp1=ftemp

!-------------------------------------------------------------
! F^n_i+1/2 (first order upwind flux)
!-------------------------------------------------------------
go to 1234
do jj=-1,np-1
 B2=B(jj,xk)
 !print*,'jj,B_jj+1/2,ftemp(jj)=',jj,B2,ftemp(jj)
 if(B2>0._8)then
  pflux(jj)=B2*ftemp(jj+1)
 else
  !temp=B2*ftemp(jj)
  !print*,'CHECKPOINT3,jj,B2,ftemp(jj),temp=',jj,B2,ftemp(jj),temp
  pflux(jj)=B2*ftemp(jj)
  !print*,'CHECKPOINT4'
 end if
 !print*,'jj,pflux(jj)=',jj,pflux(jj)
end do

!No-flux Boundary Condition
pflux(-1)=0.
pflux(np-1)=0.
1234 continue

!Update distribution function
!fupdate=ftemp1
do jj=0,np-1
a=px(jj)*dudx(xk)/3.
a=a
!ftemp(jj)=ftemp(jj)*(1.+dt*Hp(jj)) &
!           +(dt/dp(jj))*(pflux(jj)-pflux(jj-1))/ &
!           + dt*S(jj,t)
ftemp(jj)=ftemp(jj)+0.5*(dt/dp(jj))*a*( &
           (1._8-sign(1._8,a))*(ftemp1(jj+1)-ftemp1(jj)) &
          +(1._8+sign(1._8,a))*(ftemp1(jj)-ftemp1(jj-1)) )

!if(abs((dt/dp(jj))*a)>1.)then
!print*,'jj,(dt/dp(jj))*a=',jj,(dt/dp(jj))*a
!end if
!if(abs((dt/dp(jj))*a)>0.9999) STOP
!end if

end do

end subroutine momentumAdvect

end module FokkerPlanck_mod
