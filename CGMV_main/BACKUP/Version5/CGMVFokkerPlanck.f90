!Coarse Grain Momentum Volume(CGMV) algorithm for Isotropic Cosmic Ray-Fokker Planck Equation

module FokkerPlanck_mod
use constants_mod
implicit none

integer,parameter::subcycleOption=0 !0: off, 1:on
real*8::fexact(0:np-1),dt2,dtf
real*8,parameter::cour2=0.5
integer::nt2


contains

subroutine CRinit()

integer::i,j
real*8::pres,Pc1(0:nx)
real*8::q2(-1:np,-1:nx),q3(-1:np,-1:nx),qtemp

open(unit=10,file='fp.txt')
open(unit=13,file='fx.txt')

open(unit=110,file='ng_p.txt')
open(unit=111,file='ng_x.txt')


!----------------------------------------------
!Initialization
!----------------------------------------------
!Momentum space bounds
pmin=1.D-1
pmax=1.D25

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


!----------------------------------------------------------------------
!Piece-wise power law interpolation: f_(p),j= f_i-1,j*(p/p_i-1)^q_i,j
!----------------------------------------------------------------------
do i=1,np-1
 do j=0,nx-1
  !Intrabin Spatial Index
  !q(i,j)=-log(f(i-1,j)/f(i,j))/log(px(i-1)/px(i))
  q(i,j)=6.0
 end do
   q(i,-1)=q(i,0)
   q(i,nx)=q(i,nx-1)
end do 


!Initialize distribution function
do i=0,np-1
 do j=0,nx-1
  !if(px(i)>1.D-4)then
  f(i,j)=px(i)**(-q(i,j))!*exp(-((xmin+j*dx)/(10.*dx))**2. )
  !else 
  !f(i,j)=1.d-25
  !end if
  

  !Gaussian in momentum space, uniform in physical space
  !f(i,j)=1000.*exp(-( (pmin+i*dp(i)-0.5*(pmax-pmin))/(10*dp(i)) )**2. )

  !Gaussian in physicsl space, uniform in momentum space
  !f(i,j)=1000.*exp(-( (xmin+j*dx-0.5*(xmax-xmin))/(50.*dx) )**2. )
  
  !Power law in moment space, uniform in phsycial space  
  !if(i<52)then
  !  f(i,j)=50.*px(i)
  !else 
  !  f(i,j)=50.*px(51)*(px(i)/px(51))**(-5.)
  !end if

  !f(i,j)=1.D-20

 end do
 f(i,-1)=f(i,0)
 f(i,nx)=f(i,nx-1)
end do


do i=1,np-1
 do j=0,nx-1
   dw=px(i)/px(i-1)
   n(i,j)=f(i-1,j)*(px(i-1)**3.)/(3.-q(i,j))
   n(i,j)=n(i,j)*(dw**(3.-q(i,j))-1.)
   g(i,j)=f(i-1,j)*(px(i-1)**4.)/(4.-q(i,j))
   g(i,j)=g(i,j)*(dw**(4.-q(i,j))-1.)

 end do
 n(i,-1)=n(i,0)
 n(i,nx)=n(i,nx-1)
 g(i,-1)=g(i,0)
 g(i,nx)=g(i,nx-1)
end do 

!Normalize Distribution function such that CR pressure is equal to downstream fluid pressure
Pc=0.
do j=0,nx-1
 do i=1,np-1
  Pc(j)=Pc(j)+pi*(4.*cl/3.)*g(i,j)
  !Pc(j)=Pc(j)+(4.*cl/3.)*(px(i)**3.)*f(i,j)*dp2(i)
 end do
end do

!pres=(gam-1.)*(u1(nx,3)-0.5*u1(nx,2)*u1(nx,2)/u1(nx,1))

do i=0,np-1
 do j=0,nx-1
   pres=(gam-1.)*(u1(j,3)-0.5*u1(j,2)*u1(j,2)/u1(j,1))
   f(i,j)=f(i,j)*(pres/Pc(j))*100.
  end do
  f(i,-1)=f(i,0)
  f(i,nx)=f(i,nx-1)
end do

!------------------------------------------------
!CR distribution function moments: n_i,j, g_i,j
!------------------------------------------------
do i=1,np-1
 do j=0,nx-1
   dw=px(i)/px(i-1)
   n(i,j)=f(i-1,j)*(px(i-1)**3.)/(3.-q(i,j))
   n(i,j)=n(i,j)*(dw**(3.-q(i,j))-1.)

   g(i,j)=f(i-1,j)*(px(i-1)**4.)/(4.-q(i,j))
   g(i,j)=g(i,j)*(dw**(4.-q(i,j))-1.)
 end do 
 n(i,-1)=n(i,0)
 n(i,nx)=n(i,nx-1)
 g(i,-1)=g(i,0)
 g(i,nx)=g(i,nx-1)
end do 

do j=0,nx-1
  write(13,*) xmin+j*dx,f(5,j)
end do
do j=0,np-1
 write(10,*) px(j),(px(j)**4.)*f(j,1)
end do

do i=1,np-1
 write(110,*) i,px(i),n(i,1),g(i,1),q(i,1)
end do

do j=0,nx-1
 write(111,*) j,xmin+dx*j,n(10,j),g(10,j),q(10,j)
end do

!Compute CR Pressure
Pc=0.
do j=0,nx-1
  do i=1,np-1
    Pc(j)=Pc(j)+pi*(4.*cl/3.)*g(i,j)
  end do
end do
Pc(-1)=Pc(0)
Pc(nx)=Pc(nx-1)

!do j=0,nx-1
!   pres=(gam-1.)*(u1(j,3)-0.5*u1(j,2)*u1(j,2)/u1(j,1))
!  print*,'j,Pc,Pg=',j,Pc(j),pres
!end do

go to 102
!Index of Reconstructed distribution function
do j=0,nx-1   
 do i=1,np-1
  dw=px(i)/px(i-1)

  !Update Spectral Index
  call rootSolveNewton(8._8,g(i,j)/n(i,j)/px(i-1),dw,qtemp)  
  q2(i,j)=qtemp
  call rootSolveSecant(8._8,g(i,j)/n(i,j)/px(i-1),dw,qtemp)
  q3(i,j)=qtemp
 end do
end do


print*,'Spectral index before and after reconstruction:'
do i=1,np-1
 print*,'i,q, q_Newton, q_Secant=',i,q(i,250),q2(i,250),q3(i,250)
end do
102 continue

print*,'Done initializing CR distribution.'

end subroutine

subroutine CRevolve()
!Evolve Distribution Function
integer::i,j,k
real*8::Fn(-1:np,0:nx-1),Fp(-1:np,0:nx-1),r(0:nx-1),qtemp
real*8::ftemp(0:nx-1),ntemp(-1:np,-1:nx),gtemp(-1:np,-1:nx)
real*8::ntot,gtot,ftot,fluxsum1,fluxsum2,a,fnL,fnR,fgL,fgR


if(subcycleOption==1)then
!Momentum advection time step
!dt2=dt*2.
dt2=dt
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

else if(subcycleOption==0)then
 nt2=0
 dt2=dt!2._8*dt
end if

!------------------------------------------
!Momentum Space Fluxes
!------------------------------------------ 
do i=1,nt2+1
  !print*,'nt2,Subcycle step=,',nt2,i
  !print*,'dt,dt2,dtf=',dt,dt2,dtf
  
  do j=0,nx-1
   !---------------------------
   !Upwinded flux
   !---------------------------
   do k=1,np-2
    dw=px(k)/px(k-1)
    if(dudx(j)<=0._8)then
     Fn(k,j)=f(k-1,j)*dw**(-q(k,j))
     Fn(k,j)=Fn(k,j)*(1./3.)*dudx(j)*(px(k)**3.)
    else
     Fn(k,j)=f(k,j)
     Fn(k,j)=Fn(k,j)*(1./3.)*dudx(j)*(px(k)**3.)
    end if
    Fp(k,j)=px(k)*Fn(k,j)

    !Momentum Space Diffusion contribution
    !if(q(k,j)>0._8)then
    !  Fn(k,j)=Fn(k,j)-q(k,j)*px(k)*D(k)*f(k-1,j)*dw**(-q(k,j))
    !else
    !  Fn(k,j)=Fn(k,j)-q(k+1,j)*px(k)*D(k)*f(k,j)
    !end if
   end do
   !No-flux momentum boundary conditions
   Fn(0,j)=0._8
   Fn(np-1,j)=0._8
   Fp(0,j)=0._8
   Fp(np-1,j)=0._8
  end do
!------------------------------------------
!Physical Space Evolution
!------------------------------------------  

!Spatial Advection (first order upwind)
if(advectionOption==2)then

ntemp=n
gtemp=g
do j=1,np-1
 do k=0,nx-1
   a=u2(k,2)/u2(k,1)
   if(a>=0._8)then
     fnL=a*ntemp(j,k-1)
     fnR=a*ntemp(j,k)
     fgL=a*gtemp(j,k-1)
     fgR=a*gtemp(j,k)
   else
     fnL=a*ntemp(j,k)
     fnR=a*ntemp(j,k+1)
     fgL=a*gtemp(j,k)
     fgR=a*gtemp(j,k+1)
   end if 
   
  !if(k==0)then
  !  fnL=0.
  !  fgL=0.
  ! end if
  ! if(k==nx-1)then
  !  fnR=0.
  !  fgR=0.
  ! end if
   n(j,k)=ntemp(j,k)-(dt2/dx)*(fnR-fnL)
   g(j,k)=gtemp(j,k)-(dt2/dx)*(fgR-fgL)
 end do
 !Open boundaries
 n(j,-1)=n(j,0)
 g(j,-1)=g(j,0)
 n(j,nx)=n(j,nx-1)
 g(j,nx)=g(j,nx-1)
end do

end if

!print*,'CHECKPOINT1'


!Spatial Diffusion+Source Terms 	
do j=1,np-1
  ftemp=0.

  !Evolve zeroth moment: n_i,j
  do k=0,nx-1
   ftemp(k)=n(j,k)
  end do  

  call computeSpatialCoefficients(j,1)

  do k=0,nx-1
    r(k)=n(j,k)+dt2*(Fn(j,k)-Fn(j-1,k))+dt2*Qi(j,k,1)
   if(advectionOption==2)then
    r(k)=r(k)-dt2*dudx(k)*n(j,k)
   end if
  end do

  call tridiagSolveSpace(ftemp,r)
   

  do k=0,nx-1
    n(j,k)=ftemp(k)   
  end do

  !Evolve first moment: g_i,j
  do k=0,nx-1
   ftemp(k)=g(j,k)
  end do  

  call computeSpatialCoefficients(j,2)
  do k=0,nx-1
    r(k)=g(j,k)+dt2*(Fp(j,k)-Fp(j-1,k)) &
         +dt2*(-(1./3.)*dudx(k))*g(j,k) &
         +dt2*(q(j,k)*D(j)/(px(j)*px(j)))*g(j,k)+dt2*Qi(j,k,2) 
    if(advectionOption==2)then
     r(k)=r(k)+dt2*(-dudx(k))*g(j,k)  
   end if
  end do
  call tridiagSolveSpace(ftemp,r)

  do k=0,nx-1
   g(j,k)=ftemp(k)   
  end do

  !Continuous bundaries 
  n(j,-1)=n(j,0)
  g(j,-1)=g(j,0)
  n(j,nx)=n(j,nx-1)
  g(j,nx)=g(j,nx-1)
end do

 if(i==nt2)then
  dt2=dtf
 end if

end do

!print*,'CHECKPOINT2'

!----------------------------------------------
!Reconstruct distribution function from moments
!----------------------------------------------
do j=0,nx-1   
 do k=1,np-1
  dw=px(k)/px(k-1)

  !Update Spectral Index

  !call rootSolveNewton(q(k,j),g(k,j)/n(k,j)/px(k-1),dw,qtemp)  
  
  call rootSolveSecant(q(k,j),g(k,j)/n(k,j)/px(k-1),dw,qtemp)
  !print*,'CHECKPOINT2.5'
  q(k,j)=qtemp


  !Update Distribution function bin edge values
  f(k-1,j)=(3.-q(k,j))*n(k,j)/(px(k-1)**3.)
  f(k-1,j)=f(k-1,j)/(dw**(3.-q(k,j))-1.)
 end do
 dw=px(np-1)/px(np-2)
 f(np-1,j)=f(np-2,j)*dw**(-q(np-1,j)) 
end do

do j=0,np-1
  f(j,-1)=f(j,0)
  f(j,nx)=f(j,nx-1)
end do


!-----------------------------------------------
!Update CR pressure
!-----------------------------------------------
Pc=0.
do j=0,nx-1
  do k=1,np-1
    Pc(j)=Pc(j)+pi*(4.*cl/3.)*g(k,j)
  end do
end do
Pc(-1)=Pc(0)
Pc(nx)=Pc(nx-1)

end subroutine CRevolve


subroutine computeSpatialCoefficients(pk,option)

!Input variables
integer::pk,xk,option
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

!-------------------------------------
!No-flux/Continuous Boundary condition 
!-------------------------------------

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
  Atx(jj)=dt2/dx
  Atx(jj)=Atx(jj)*(Wminus(jj-1)-Cp(jj-1))
  if(methodType==2)then
   Atx(jj)=Atx(jj)*0.5 
  end if
  
  Ctx(jj)=dt2/dx
  Ctx(jj)=-Ctx(jj)*(Wplus(jj)+Cp(jj))
  if(methodType==2)then
   Ctx(jj)=Ctx(jj)*0.5
  end if

   Btx(jj)= Wminus(jj)-Wplus(jj-1)-Cp(jj)-Cp(jj-1)
   Btx(jj)=Btx(jj)*dt2/dx
   Btx(jj)=-Btx(jj)+1.
  if(methodType==2)then 
   Btx(jj)=0.5*(Btx(jj)+1.)
  end if

  !print*,'jj,Atx,Btx,Ctx=',jj,Atx(jj),Btx(jj),Ctx(jj)
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

subroutine tridiagSolveSpace(ftemp,r) !from Numerical Recipes book

!Input variables
real*8::ftemp(0:nx-1),r(0:nx-1)

real*8::bet,gam(0:nx-1)
integer::kk


if(methodType==2)then
 do kk=0,nx-1
   r(kk)=r(kk)+ftemp(kk)-Atx(kk)*ftemp(kk-1)-Btx(kk)*ftemp(kk)-Ctx(kk)*ftemp(kk+1)
 end do
end if

bet=Btx(0)
ftemp(0)=r(0)/bet

!Upper triangular decomposition and forward substitution
do kk=1,nx-1
  gam(kk)=Ctx(kk-1)/bet
  bet=Btx(kk)-Atx(kk)*gam(kk)
  ftemp(kk)=(r(kk)-Atx(kk)*ftemp(kk-1))/bet
end do

!Back substitution
do kk=nx-2,0,-1
 ftemp(kk)=ftemp(kk)-gam(kk+1)*ftemp(kk+1)
end do

end subroutine tridiagSolveSpace

!------------------------------------------------------------------!
!------------------------------------------------------------------!
!Momentum Diffusion co-efficient
function D(pk) result(fx)

integer,intent(in)::pk
real*8::fx,Dpp

 Dpp=0._8
 fx=Dpp

end function D

!Spatial Advection co-efficient
function E(xk) result(ex)
integer,intent(in)::xk
real*8::ex


if(advectionOption==1)then
  ex=-u2(xk,2)/u2(xk,1)
else if(advectionOption==2)then
  ex=0._8
end if

end function E

!Spatial Diffusion co-efficient
function Kxx(xk,pk) result(kx)

integer,intent(in)::xk,pk
real*8::kx

!kx=0.0005
kx=0.05/u2(xk,1)
!kx=0.1*(px(pk)**0.51)
end function Kxx

!CR particle number Injection term
function Qi(pk,xk,option) result(sx)
 
integer,intent(in)::pk,xk,option
real*8::sx

real*8,parameter::alpha=2.
real*8::eps=0.001
real*8::mp=0.1

real*8::cs2,rho1,rho2,pinj,us,pres,wx
real*8::ut1,ut2
integer::pin

!------------------------------
!Flux fraction injection model
!------------------------------
!Shock Speed
!ut1=-2.*sqrt(gam)    !u1(shockCell(1),2)/u1(shockCell(1),1)
!ut2=0.               !u1(shockCell(1)-2,2)/u1(shockCell(1)-2,1)
ut1=u2(shockCell(1),2)/u2(shockCell(1),1)
ut2=u2(shockCell(1)-2,2)/u2(shockCell(1)-2,1)
us=(1./3.)*(4.*ut2-ut1)   !assuming strong shock
!pre-shock density
!rho1=1.!u1(shockCell(1),1)
rho1=u2(shockCell(1),1)
!Post-shock density
rho2=u2(shockCell(1)-1,1)
!Post-shock sound speed
pres=(gam-1.)*(u2(shockCell(1)-1,3) &
-0.5*u2(shockCell(1)-1,2)*u2(shockCell(1)-1,2)/u2(shockCell(1)-1,1))

cs2=sqrt(gam*pres/rho2)

!Compute injection momentum
pinj=1.5!alpha*cs2
pin=np*log(pinj/pmin)/log(pmax/pmin)

!print*,'pinj,pin=',pinj

!Spatial weight function
!wx=(dx*xk-dx*(shockCell(1)-1))/(10.*dx)
!wx=exp(-wx**2.)/sqrt(3.141592*(10.*dx)**2.)
if(xk==shockCell(1)-1)then
  wx=1.
else
  wx=0.
end if

!Injection source term
if(pinj<= px(pk) .and. pinj>= px(pk-1))then
 sx=wx*eps*rho1*us/mp 
else
 sx=0._8
end if


if(option==2)then
sx=sx*pinj
end if

sx=0.

end function Qi

function fluidS(xk) result(fx)
integer,intent(in)::xk
real*8::fx

real*8,parameter::alpha=2.
real*8::eps=0.001
real*8::mp=0.1

real*8::cs2,rho1,rho2,pinj,us,pres,wx
real*8::ut1,ut2
integer::pin

!------------------------------
!Flux fraction injection model
!------------------------------
!Shock Speed
ut1=u2(shockCell(1),2)/u2(shockCell(1),1)
ut2=u2(shockCell(1)-2,2)/u2(shockCell(1)-2,1)
us=(1./3.)*(4.*ut2-ut1)     !(1./3.)*2.*sqrt(gam)
!pre-shock density
rho1=1.!u1(shockCell(1),1)
!Post-shock sound speed
!pres=55.78043
!pres=(gam-1.)*(u1(shockCell(1)-2,3) &
!-0.5*u1(shockCell(1)-2,2)*u1(shockCell(1)-2,2)/u1(shockCell(1)-2,1))
!rho2=u1(shockCell(1)-2,1)
pres=(gam-1.)*(u2(shockCell(1)-1,3) &
-0.5*u2(shockCell(1)-1,2)*u2(shockCell(1)-1,2)/u2(shockCell(1)-1,1))
rho2=u2(shockCell(1)-1,1)
cs2=sqrt(gam*pres/rho2)

!Compute injection momentum
!pinj=lambda*sqrt(gam*pres/ 3.976468)
pinj=1.5!alpha*cs2
pin=np*log(pinj/pmin)/log(pmax/pmin)

!Spatial weight function
!wx=(dx*xk-dx*(shockCell(1)-1))/(10.*dx)
!wx=exp(-wx**2.)/sqrt(3.141592*(10.*dx)**2.)
if(xk==shockCell(1)-1)then
  wx=1.
else
  wx=0.
end if

fx=0.5*eps*wx*(pinj**2.)*rho1*us

end function fluidS
!------------------------------------------------------------------!
!------------------------------------------------------------------!

subroutine rootSolveNewton(qi0,gnp,dw,qi)
real*8::qi0,gnp,dw,qi
integer::niter=20
real*8::q0
real*8::itertol=1.D-6
integer::i,iterCount=0

!Intial guess
q0=qi0

do i=1,niter
  qi=q0-psi2(q0,gnp,dw)/psip2(q0,gnp,dw)
  if(debug==1)then
   print*,'iteration #',i+iterCount*20
   print*,'q0,qi=',q0,qi
  end if 

  if(abs(qi-q0)<itertol)then
    if(debug==1)then
    print*,'Iterations converged.'
    print*,'qi0,qi, #iteration=',qi0,qi,i+iterCount*20
    end if 
    EXIT
  end if
  if(abs(qi-q0)>itertol .and. i==niter-1 .and. iterCount<5)then
    niter=niter+20
    iterCount=iterCount+1
  end if

  if(abs(qi-q0)>itertol .and. i==niter-1 .and. iterCount==5)then
     print*,'Iterations fail to converge...terminating suboutine.'
     qi=qi0
     EXIT
  end if

  q0=qi
end do


end subroutine rootSolveNewton

subroutine rootSolveSecant(qi0,gnp,dw,qi)
real*8::qi0,gnp,dw,qi
integer::niter=20
real*8::q0,q1
real*8::itertol=1.D-6
integer::i,iterCount=0

!Intial guesses
q0=qi0
q1=1.00001*qi0

do i=1,niter
  qi=q1-psi2(q1,gnp,dw)/((psi2(q1,gnp,dw)-psi2(q0,gnp,dw))/(q1-q0))
  if(debug==1)then
   print*,'iteration #',i+iterCount*20
   print*,'q0,q1,qi=',q0,q1,qi
  end if 

  if(abs(qi-q1)<itertol)then
    if(debug==1)then
    print*,'Iterations converged.'
    print*,'qi0,qi, #iteration=',qi0,qi,i+iterCount*20
    end if 
    EXIT
  end if
  if(abs(qi-q1)>itertol .and. i==niter-1 .and. iterCount<5)then
    niter=niter+20
    iterCount=iterCount+1
  end if

  if(abs(qi-q1)>itertol .and. i==niter-1 .and. iterCount==10)then
     print*,'Iterations fail to converge...terminating suboutine.'
     qi=qi0
     EXIT
  end if

  q0=q1
  q1=qi

end do

end subroutine rootSolveSecant


function psi2(qi,gnp,dw) result(fx)
real*8,intent(in)::qi,gnp,dw
real*8::fx,t1,t2,t3,t4

t1=3._8-qi
t2=(dw**(4._8-qi))-1.
t3=4._8-qi
t4=(dw**(3._8-qi))-1.

fx=(t1*t2/(t3*t4))-gnp

end function psi2

function psip2(qi,gnp,dw) result(fx)
real*8,intent(in)::qi,gnp,dw
real*8::fx,t1,t2,t3,t4

t1=1./(4._8-qi)/((dw**(3._8-qi))-1.)
t2=psi2(qi,gnp,dw)+gnp
t3=((dw**(3._8-qi))-1.)+(4._8-qi)*log(dw)*(dw**(3._8-qi))
t4=((dw**(4._8-qi))-1.)+(3._8-qi)*log(dw)*(dw**(4._8-qi))

fx=t1*(t2*t3-t4)

end function psip2


end module FokkerPlanck_mod
