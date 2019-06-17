!Chang-Cooper(1970) algorithm for Isotropic Cosmic Ray-Fokker Planck Equation

!Note1: Momentum (p) is in units of m_p c (where m_p is the CR particle mass 
!and c is the speed of light)

!Note2: For DSA modified shock simulation, require dx << min(D_diff) ~ min(K_diff)/u_sh 
!to get accurate convergence of solution. (D_diff is the diffusion length)

!Note3: Issue with sub-cycling 

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
pmin=1.D-2
pmax=1.D40

print*,'np,pmin,pmax=',np,pmin,pmax

nparticle=0.

!Generate momentum mesh points
do i=-2,np+3
  if(momentumMeshType==1)then
   px(i)=pmin+i*(pmax-pmin)/np
  else if(momentumMeshType==2)then
   px(i)=pmin*(pmax/pmin)**(real(i)/real(np-1))
  end if
  !print*,'i,px(i)=',i,px(i)
end do


!Momentum mesh spacing
do i=-1,np+2
 dp(i)=0.5*(px(i+1)-px(i-1))
 dp2(i)=px(i)-px(i-1)
 !print*,'i,p,dp=',i,px(i),dp(i)
end do

!Check diffusion length and spatial resolution requirement
xd=Kxx(0,0)
do i=0,np-1
 do j=0,nx-1
  xd=min(xd,Kxx(j,i))
 end do
end do

print*,'xd,dx',xd,dx!,dx/xd

do i=0,np-1
 do j=0,nx-1
  !Initial profile
  !if(px(i)>=10.)then
  ! f(i,j)=px(i)**(-5.)!*exp(-((xmin+j*dx)/(10.*dx))**2. )
  !end if
  f(i,j)=px(i)**(-4.5)
  !f(i,j)=f(i,j)*exp(-( (xmin+j*dx-0.5*(xmax-xmin))/(50.*dx) )**2. )
  !f(i,j)=0._8
  !f(i,j)=1000.*exp(-( (xmin+j*dx-0.5*(xmax-xmin))/(5.*dx) )**2. )!
  !f(i,j)=1000.*exp(-( (pmin+i*dp(i)-0.5*(pmax-pmin))/(20.*dp(i)) )**2. ) 
 end do
 f(-1,j)=f(0,j)  
 f(np,j)=f(np-1,j)
end do



!Normalize Distribution function such that CR pressure is equal to downstream fluid pressure
do j=0,nx-1
 do i=1,np-1
  Pc(j)=Pc(j)+(4.*cl/3.)*(px(i)**3.)*f(i,j)*dp2(i)
 end do
end do

pres=(gam-1.)*(u1(nx,3)-0.5*u1(nx,2)*u1(nx,2)/u1(nx,1))

do i=0,np-1
 do j=0,nx-1
  f(i,j)=f(i,j)*(pres/Pc(j))!*100.
 end do
end do


fp=0.
do i=0,np-1
 do j=0,nx-1
  fp(i)=fp(i)+dx*f(i,j)
 end do
 write(10,*) px(i),f(i,0)
 !print*,'p,f(p)',px(i),fp(i)
end do

fx=0.
do j=0,nx-1
  do i=0,np-1
  fx(j)=fx(j)+dp(i)*f(i,j)
 end do
 write(13,*) xmin+j*dx,f(25,j)
 !print*,'x,f(x)',xmin+j*dx,fx(j)
end do


!Compute CR Pressure
Pc=0.
do j=0,nx-1
  do i=1,np-1
    Pc(j)=Pc(j)+(4.*cl/3.)*(px(i)**3.)*f(i,j)*dp2(i)
  end do
end do
Pc(-1)=Pc(0)
Pc(-2)=Pc(0)
Pc(nx)=Pc(nx-1)
Pc(nx+1)=Pc(nx-1)


!do j=0,nx-1
!   pres=(gam-1.)*(u1(j,3)-0.5*u1(j,2)*u1(j,2)/u1(j,1))
!  print*,'j,Pc,Pg=',j,Pc(j),pres
!end do


print*,'Done initializing CR distribution.'

end subroutine

subroutine CRevolve()
!Evolve Distribution Function
integer::i,j,k
real*8::pres,f2(-1:np,-1:nx)

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
  dtf=0.
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

  call computeSpatialCoefficients(j)
  call tridiagSolveSpace(nx,ftemp2)

  ftemp2(-1)=ftemp2(0)  
  ftemp2(nx)=ftemp2(nx-1)

  do k=-1,nx
  f(j,k)=ftemp2(k)   
  end do
 end do

 if(i==nt2) dt3=dtf
end do


!-----------------------------------------------
!Update CR pressure
!-----------------------------------------------
Pc=0.
do j=0,nx-1
  do k=1,np-1
    Pc(j)=Pc(j)+(4.*cl/3.)*(px(k)**3.)*f(k,j)*dp2(k)
  end do
end do
Pc(-1)=Pc(0)
Pc(-2)=Pc(0)
Pc(nx)=Pc(nx-1)
Pc(nx+1)=Pc(nx-1)

!do j=-1,nx
!   pres=(gam-1.)*(u1(j,3)-0.5*u1(j,2)*u1(j,2)/u1(j,1))
!  print*,'j,Pc,Pg=',j,Pc(j),pres
!end do


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

!for no-flux boundary on the left
Wminus(-1)=0.
Wplus(-1)=0.
Cp(-1)=0.

!for no-flux boundary on the right
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
end do

!for open boundary on the left 
!Atx(0)=0.
!Btx(0)=1.
!Ctx(0)=0.

!for open boundary on the right 
Atx(nx-1)=0.
Btx(nx-1)=1.
Ctx(nx-1)=0.

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

ex=-u2(xk,2)/u2(xk,1)

end function E

!Spatial Diffusion co-efficient
function Kxx(xk,pk) result(kx)

integer,intent(in)::xk,pk
real*8::kx

!kx=0.05
kx=0.05/u2(xk,1)
!kx=10.*(px(pk)**0.51)


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
ut1=-10.*sqrt(gam)    !u1(shockCell(1),2)/u1(shockCell(1),1)
ut2=0.               !u1(shockCell(1)-2,2)/u1(shockCell(1)-2,1)
us=(1./3.)*(4.*ut2-ut1)     !(1./3.)*2.*sqrt(gam)
!pre-shock density
rho1=1.!u1(shockCell(1),1)
!Post-shock sound speed
!pres=55.78043
!pres=(gam-1.)*(u1(shockCell(1)-2,3) &
!-0.5*u1(shockCell(1)-2,2)*u1(shockCell(1)-2,2)/u1(shockCell(1)-2,1))
!rho2=u1(shockCell(1)-2,1)
pres=(gam-1.)*(u1(shockCell(1)-1,3) &
-0.5*u1(shockCell(1)-1,2)*u1(shockCell(1)-1,2)/u1(shockCell(1)-1,1))
rho2=u1(shockCell(1)-1,1)
cs2=sqrt(gam*pres/rho2)

!Compute injection momentum
!pinj=lambda*sqrt(gam*pres/ 3.976468)
pinj=alpha*cs2
pin=np*log(pinj/pmin)/log(pmax/pmin)

!print*,'pinj,pin=',pinj,pin


!sx=0.
!injection at shock location
!if(xk==shockCell(1)-1)then
!  sx=(px(pk)-100.)/(10.)
!  sx=exp(-sx**2.)
if(pk==pin .and. pin>0 .and. pin <np-1 )then
 
 !Spatial weight function
 wx=(dx*xk-dx*(shockCell(1)-1))/(4.*dx)
 wx=exp(-wx**2.)/sqrt(3.141592*(4.*dx)**2.)
 !if(xk==shockCell(1)-1)then
 !  wx=1.
 !else
 !  wx=0.
 !end if
 !Injection source term
 sx=wx*eps*rho1*us/(4.*3.141592*px(pk)*px(pk))/mp 
else
  sx=0._8
end if

sx=0._8

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
real*8::fx

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
ut1=-10.*sqrt(gam)    !u1(shockCell(1),2)/u1(shockCell(1),1)
ut2=0.               !u1(shockCell(1)-2,2)/u1(shockCell(1)-2,1)
us=(1./3.)*(4.*ut2-ut1)     !(1./3.)*2.*sqrt(gam)
!pre-shock density
rho1=1.!u1(shockCell(1),1)
!Post-shock sound speed
!pres=55.78043
!pres=(gam-1.)*(u1(shockCell(1)-2,3) &
!-0.5*u1(shockCell(1)-2,2)*u1(shockCell(1)-2,2)/u1(shockCell(1)-2,1))
!rho2=u1(shockCell(1)-2,1)
pres=(gam-1.)*(u1(shockCell(1)-1,3) &
-0.5*u1(shockCell(1)-1,2)*u1(shockCell(1)-1,2)/u1(shockCell(1)-1,1))
rho2=u1(shockCell(1)-1,1)
cs2=sqrt(gam*pres/rho2)

!Compute injection momentum
!pinj=lambda*sqrt(gam*pres/ 3.976468)
pinj=alpha*cs2
pin=np*log(pinj/pmin)/log(pmax/pmin)

!Spatial weight function
wx=(dx*xk-dx*(shockCell(1)-1))/(4.*dx)
wx=exp(-wx**2.)/sqrt(3.141592*(4.*dx)**2.)
!if(xk==shockCell(1)-1)then
!  wx=1.
!else
!  wx=0.
!end if

fx=0.5*eps*wx*(pin**2.)*rho1*us

fx=0._8

end function fluidS


end module FokkerPlanck_mod
