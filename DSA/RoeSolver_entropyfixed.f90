!Roe's Approximate (Linear) Riemann Solver (muscl-hancock 2ND ORDER)
module RoeSolver_entropyfixed_mod
use constants_mod
implicit none

real*8::rhoL,uL,pL,aL,HL,FL(3) !left state
real*8::rhoR,uR,pR,aR,HR,FR(3) !right state 
real*8::lambda1L,lambda1R,lambda3L,lambda3R
real*8::pS,rhoSL,rhoSR,uS,aSL,aSR
real*8::smax,smax1
real*8::sL,sR

contains 

subroutine init()

integer::j,k
real*8::rhoL,uL,pL !left state
real*8::rhoR,uR,pR!right state 
!------------------------------------------------------
!initialization
!------------------------------------------------------

!open(unit=10,file='input.txt')
open(unit=11,file='output_fluid.txt')
open(unit=12,file='dt.txt')

!Discountinuity in initial condition at the middle of [xmin,xmax]

dx=(xmax-xmin)/nx

gam=5./3.
gamil=1./(gam+1.)
gamul=1./(gam-1.)
gamel=(gam-1.)/(gam+1.)
gamee=(gam-1.)/(2.*gam)
gamuu=(gam+1.)/(2.*gam)

smax=0.

!read(unit=10,fmt=*) rhoL,uL,pL 
!read(unit=10,fmt=*) rhoR,uR,pR

!rhoL=1.
!rhoR=rhoL
!pL=1.
!pR=pL
!uL=-10.*sqrt(gam*pL/rhoL)
!uR=uL

rhoL=1.
rhoR=rhoL
pL=0.001
pR=pL
uL=-40.*sqrt(gam*pL/rhoL)
uR=uL

!rhoL=1.
!rhoR=rhoL
!uL=1.
!uR=uL
!pL=(1./gam)*(1./30.)**2.!0.1
!pR=pL


do j=-2,nx+1
  if(j<0) then
    u1(j,1)=rhoL
    u1(j,2)=rhoL*uL
    u1(j,3)=0.5*rhoL*uL*uL+pL*gamul
  else
    u1(j,1)=rhoR
    u1(j,2)=rhoR*uR
    u1(j,3)=0.5*rhoR*uR*uR+pR*gamul
  end if
end do

u2=u1

print*,'Left State (rho,u,p) =',rhoL,uL,pL
print*,'Right State (rho,u,p) =',rhoR,uR,pR


!close(unit=10)

end subroutine init

subroutine computeFlux()	

real*8::delrho,delp,delu
real*8::dens,vel,pres
real*8::rhot,ut,Ht,at,K1(3),K2(3),K3(3),alpha(3)
real*8::qL(-2:nx+1,3),qR(-2:nx+1,3),slope(-2:nx+1,3)
real*8::ubL(-2:nx+1,3),ubR(-2:nx+1,3)


integer::j,k

!--------------------------------------------------
!Limited Slopes: delta_j
!--------------------------------------------------
do j=-1,nx
  do k=1,3
   slope(j,k)=limiter(u1(j,k)-u1(j-1,k),u1(j+1,k)-u1(j,k))
  end do
end do

!----------------------------------------------
!Boundary Extrapolated Values: uL_i, uR_i
!----------------------------------------------
do j=-1,nx
  do k=1,3
   qL(j,k)=u1(j,k)-0.5*slope(j,k)
   qR(j,k)=u1(j,k)+0.5*slope(j,k)
  end do
end do

!-------------------------------------------------------------------
!Evolve Boudary Extrapolated values by a half-step: ubarL_i, ubarR_i
!-------------------------------------------------------------------
do j=-1,nx
  FL(1)=qL(j,2)
  FL(2)=(gam-1.)*qL(j,3)+0.5*(3.-gam)*(qL(j,2)**2.)/qL(j,1)
  FL(3)=gam*(qL(j,2)*qL(j,3)/qL(j,1))-0.5*(gam-1.)*(qL(j,2)**3.)/(qL(j,1)**2.)

  FR(1)=qR(j,2)
  FR(2)=(gam-1.)*qR(j,3)+0.5*(3.-gam)*(qR(j,2)**2.)/qR(j,1)
  FR(3)=gam*(qR(j,2)*qR(j,3)/qR(j,1))-0.5*(gam-1.)*(qR(j,2)**3.)/(qR(j,1)**2.)

  do k=1,3
    ubL(j,k)=qL(j,k)+0.5*(dt/dx)*(FL(k)-FR(k))
    ubR(j,k)=qR(j,k)+0.5*(dt/dx)*(FL(k)-FR(k))
  end do
end do


do j=-1,nx-1
!---------------------------------
!w_j+1/2(0)
!---------------------------------

go to 112
!set left and right states for local Riemann problem at j+1/2 cell interface
rhoL=max(densmin,ubR(j,1))
uL=ubR(j,2)/ubR(j,1)
pL=max(premin,(gam-1.)*(ubR(j,3)-0.5*ubR(j,2)*ubR(j,2)/ubR(j,1)))
aL=sqrt(gam*pL/rhoL)
HL=0.5*uL*uL+gamul*aL*aL


rhoR=max(densmin,ubL(j+1,1))
uR=ubL(j+1,2)/ubL(j+1,1)
pR=max(premin,(gam-1.)*(ubL(j+1,3)-0.5*ubL(j+1,2)*ubL(j+1,2)/ubL(j+1,1)))
aR=sqrt(gam*pR/rhoR)
HR=0.5*uR*uR+gamul*aR*aR
112 continue

!go to 122 
rhoL=max(densmin,u1(j,1))
uL=u1(j,2)/u1(j,1)
pL=max(premin,(gam-1.)*(u1(j,3)-0.5*u1(j,2)*u1(j,2)/u1(j,1)))
aL=sqrt(gam*pL/rhoL)
HL=0.5*uL*uL+gamul*aL*aL

rhoR=max(densmin,u1(j+1,1))
uR=u1(j+1,2)/u1(j+1,1)
pR=max(premin,(gam-1.)*(u1(j+1,3)-0.5*u1(j+1,2)*u1(j+1,2)/u1(j+1,1)))
aR=sqrt(gam*pR/rhoR)
HR=0.5*uR*uR+gamul*aR*aR
!122 continue

delrho=rhoR-rhoL
delu=uR-uL
delp=pR-pL

!Modify Fluxes to include Cosmic Ray pressure
FL(1)=rhoL*uL
FL(2)=rhoL*uL*uL+pL!+Pc(j)
FL(3)=0.5*rhoL*uL*uL*uL+gam*gamul*pL*uL!+Pc(j)*uL

FR(1)=rhoR*uR
FR(2)=rhoR*uR*uR+pR!+Pc(j+1)
FR(3)=0.5*rhoR*uR*uR*uR+gam*gamul*pR*uR!+Pc(j+1)*uR

if(debug==1)then
  print*,'Cell: ',j
  print*,'Left State=',rhoL,uL,pL
  print*,'Right State=',rhoR,uR,pR
end if

!Compute star region state
call starRegion()

!-------------------------------------------
!Wave Speeds for Harten-Hyman Entropy Fix
!-------------------------------------------
lambda1L=uL-aL
lambda1R=uS-aSL

lambda3L=uS+aSR
lambda3R=uR+aR

!-------------------------------------------
!Roe-Averaged State Variables: ^rho,^u,^H,^u
!-------------------------------------------
rhot=sqrt(rhoL*rhoR)
ut=(sqrt(rhoL)*uL+sqrt(rhoR)*uR)/(sqrt(rhoL)+sqrt(rhoR))
Ht=(sqrt(rhoL)*HL+sqrt(rhoR)*HR)/(sqrt(rhoL)+sqrt(rhoR))
at=sqrt((gam-1.)*(Ht-0.5*ut*ut))

!---------------------------------------------
!Roe Jacobian Matrix Eigenvalues:^lambda_1,2,3
!---------------------------------------------
lambda(j,1)=ut-at
lambda(j,2)=ut
lambda(j,3)=ut+at

if(debug==1)then
  print*,'Roe average Eigenvalues=',lambda(j,1),lambda(j,2),lambda(j,3)
end if
!---------------------------------------------
!Roe Jacobian Matrix Eigenvectors:^K_1,2,3
!---------------------------------------------
K1(1)=1._8
K1(2)=ut-at
K1(3)=Ht-ut*at

K2(1)=1
K2(2)=ut
K2(3)=0.5*ut*ut

K3(1)=1._8
K3(2)=ut+at
K3(3)=Ht+ut*at

!---------------------------------------------
!Roe Averaged Wave Strengths:^alpha_1,2,3
!---------------------------------------------
alpha(1)=(0.5/(at*at))*(delp-rhot*at*delu)
alpha(2)=delrho-delp/(at*at)
alpha(3)=(0.5/(at*at))*(delp+rhot*at*delu)

!---------------------------------------
!Max Wave Speed
!---------------------------------------
smax=max(smax,abs(lambda(j,1)),abs(lambda(j,2)),abs(lambda(j,3)))

if(debug==1)then
  print*,'Max wave speed=',smax
end if
!-------------------------------------------------------------
!Numerical Flux: F_j+1/2
!-------------------------------------------------------------
!Check for transonic rarefaction waves
if(lambda1L<0. .and. lambda1R>0.)then
  do k=1,3
    flux(j,k)=FL(k)+(lambda1L*(lambda1R-lambda(j,1))/(lambda1R-lambda1L))&
            *alpha(1)*K1(k)  
  end do
else if(lambda3L<0. .and. lambda3R>0.)then
  do k=1,3
    flux(j,k)=FR(k)-(lambda3R*(lambda(j,3)-lambda3L)/(lambda3R-lambda3L))&
            *alpha(3)*K3(k)
  end do
else
  do k=1,3
    flux(j,k)=0.5*(FL(k)+FR(k))-0.5*(alpha(1)*abs(lambda(j,1))*K1(k)&
          +alpha(2)*abs(lambda(j,2))*K2(k)+alpha(3)*abs(lambda(j,3))*K3(k)) 
  end do
end if

!-----------------------------------
!Roe solution on t=0 line: u(x/t=0)
!-----------------------------------
do k=1,3
  u3(j,k)=0.5*(u1(j,k)+u1(j+1,k))-0.5*(alpha(1)*sign(1._8,lambda(j,1))*K1(k) &
         +alpha(2)*sign(1._8,lambda(j,2))*K2(k)+alpha(3)*sign(1._8,lambda(j,3))*K3(k)) 
end do

!ustar(j)=(u1(j,2)+alpha(3)*K3(2))/(u1(j,1)+alpha(3)*K3(1))

end do

end subroutine computeFlux


subroutine computeFluxHLL()

real*8::delrho,delp,delu
real*8::dens,vel,pres
real*8::rhot,ut,Ht,at,K1(3),K2(3),K3(3),alpha(3)
real*8::sL,sR
integer::j,k

do j=-1,nx-1
!---------------------------------
!w_j+1/2(0)
!---------------------------------
rhoL=u1(j,1)
uL=u1(j,2)/u1(j,1)
pL=(gam-1.)*(u1(j,3)-0.5*u1(j,2)*u1(j,2)/u1(j,1))

rhoR=u1(j+1,1)
uR=u1(j+1,2)/u1(j+1,1)
pR=(gam-1.)*(u1(j+1,3)-0.5*u1(j+1,2)*u1(j+1,2)/u1(j+1,1))

if(debug==1)then
  print*,'Cell: ',j
  print*,'Left State=',rhoL,uL,pL
  print*,'Right State=',rhoR,uR,pR
end if

aL=sqrt(gam*pL/rhoL)
HL=0.5*uL*uL+gamul*aL*aL
aR=sqrt(gam*pR/rhoR)
HR=0.5*uR*uR+gamul*aR*aR

delrho=rhoR-rhoL
delu=uR-uL
delp=pR-pL

!Modify Fluxes to include Cosmic Ray pressure
FL(1)=rhoL*uL
FL(2)=rhoL*uL*uL+pL
FL(3)=0.5*rhoL*uL*uL*uL+gam*gamul*pL*uL

FR(1)=rhoR*uR
FR(2)=rhoR*uR*uR+pR
FR(3)=0.5*rhoR*uR*uR*uR+gam*gamul*pR*uR



!Compute star region state
!call starRegion()

!-------------------------------------------
!Roe-Averaged State Variables: ^rho,^u,^H,^u
!-------------------------------------------
rhot=sqrt(rhoL*rhoR)
ut=(sqrt(rhoL)*uL+sqrt(rhoR)*uR)/(sqrt(rhoL)+sqrt(rhoR))
Ht=(sqrt(rhoL)*HL+sqrt(rhoR)*HR)/(sqrt(rhoL)+sqrt(rhoR))
at=sqrt((gam-1.)*(Ht-0.5*ut*ut))

!---------------------------------------------
!Roe Jacobian Matrix Eigenvalues:^lambda_1,2,3
!---------------------------------------------
lambda(j,1)=ut-at
lambda(j,2)=ut
lambda(j,3)=ut+at

if(debug==1)then
  print*,'Roe average Eigenvalues=',lambda(j,1),lambda(j,2),lambda(j,3)
end if
!---------------------------------------------
!Roe Jacobian Matrix Eigenvectors:^K_1,2,3
!---------------------------------------------
K1(1)=1._8
K1(2)=ut-at
K1(3)=Ht-ut*at

K2(1)=1
K2(2)=ut
K2(3)=0.5*ut*ut

K3(1)=1._8
K3(2)=ut+at
K3(3)=Ht+ut*at

!---------------------------------------------
!Roe Averaged Wave Strengths:^alpha_1,2,3
!---------------------------------------------
alpha(1)=(0.5/(at*at))*(delp-rhot*at*delu)
alpha(2)=delrho-delp/(at*at)
alpha(3)=(0.5/(at*at))*(delp+rhot*at*delu)

!---------------------------------------
!Max Wave Speed
!---------------------------------------
smax=max(smax,abs(lambda(j,1)),abs(lambda(j,2)),abs(lambda(j,3)))

if(debug==1)then
  print*,'Max wave speed=',smax
end if
!-------------------------------------------------------------
!Numerical Flux: F_j+1/2
!-------------------------------------------------------------

!Extremal Wave Speeds
sL=min(uL-aL,ut-at,-0.0000000001)
sR=max(uR+aR,ut+at,0.0000000001)

if(debug==1)then
  print*,'HLL solver Extremal wave speeds=',sL,sR
end if

do k=1,3
 if(sL>=0._8)then
   flux(j,k)=FL(k)
  else if(sL<0._8 .and. sR>0._8)then
   flux(j,k)=sR*FL(k)-sL*FR(k)+sL*sR*(u1(j+1,k)-u1(j,k))
   flux(j,k)=flux(j,k)/(sR-sL)
  else
   flux(j,k)=FR(k)
 end if
end do


!-----------------------------------
!Roe solution on t=0 line: u(x/t=0)
!-----------------------------------
do k=1,3
  u3(j,k)=0.5*(u1(j,k)+u1(j+1,k))-0.5*(alpha(1)*sign(1._8,lambda(j,1))*K1(k) &
         +alpha(2)*sign(1._8,lambda(j,2))*K2(k)+alpha(3)*sign(1._8,lambda(j,3))*K3(k)) 
end do


!ustar(j)=(u1(j,2)+alpha(3)*K3(2))/(u1(j,1)+alpha(3)*K3(1))

end do
end subroutine computeFluxHLL


subroutine computeWaveSpeeds()

!Local Variables
real*8::rhoB,aB,pMin,pMax,p0

!sound speeds for left and right states
aL=sqrt(gam*pL/rhoL)
aR=sqrt(gam*pR/rhoR)

!------------------------------------------------------
!compute pressure and velocity in star region
!------------------------------------------------------
rhoB=0.5*(rhoL+rhoR)
aB=0.5*(aL+aR)

pMin=min(pL,pR)
pMax=max(pL,pR)

pS=0.5*(pL+pR)+0.5*(uL-uR)*rhoB*aB

if(pS<pMax .and. pS>pMin .and. (pMax/pMin)<Q_user)then

!PV solver(solution of linearized primitive variable Euler eqns)

else if(pS<=pMin) then

!Two-Rarefaction Riemann Solver
!p_star
pS=aL+aR-0.5*(gam-1.)*(uR-uL)
pS=pS/(aL/(pL**gamee)+aR/(pR**gamee))
pS=pS**(1./gamee)
else if(pS>pMin) then

!Two-Shock Riemann Solver
!p_star
p0=max(0._8,pS)
pS=gL(p0)*pL+gR(p0)*pR-(uR-uL)
pS=pS/(gL(p0)+gR(p0))
end if

!Compute Wave Speeds

!------------------------------------------------------
!Two-Shocks
!------------------------------------------------------
if(pS>pR .and. pS>pL)then

!HLL Wave speeds
sL=uL-aL*sqrt(gamuu*(pS/pL)+gamee)
sR=uR+aR*sqrt(gamuu*(pS/pR)+gamee)
!------------------------------------------------------
!Two-Rarefactions
!------------------------------------------------------
else if(pS<=pL .and. pS<=pR)then

!HLL Wave speeds
sL=uL-aL
sR=uR+aR
!------------------------------------------------------
!Left Shock, Right Rarefaction
!------------------------------------------------------
else if(pS>pL .and. pS<=pR)then

!HLL Wave speeds
sL=uL-aL*sqrt(gamuu*(pS/pL)+gamee)
sR=uR+aR
!------------------------------------------------------
!Left Rarefaction, Right Shock
!------------------------------------------------------
else if(pS<=pL .and. pS>pR)then

!HLL Wave speeds
sL=uL-aL
sR=uR+aR*sqrt(gamuu*(pS/pR)+gamee)
end if

smax1=abs(sR)

end subroutine computeWaveSpeeds



subroutine starRegion()

!Local Variables
real*8::rhoB,aB,pMin,pMax,p0

rhoB=0.5*(rhoL+rhoR)
aB=0.5*(aL+aR)

pMin=min(pL,pR)
pMax=max(pL,pR)

pS=0.5*(pL+pR)+0.5*(uL-uR)*rhoB*aB

if(pS<pMax .and. pS>pMin .and. (pMax/pMin)<Q_user)then

!PV solver(solution of linearized primitive variable Euler eqns)
!u_star
uS=0.5*(uL+uR)+0.5*(pL-pR)/(rhoB*aB)
!rho_star
rhoSL=rhoL+(uL-uS)*(rhoB/aB)
rhoSR=rhoR+(uS-uR)*(rhoB/aB)

else if(pS<=pMin) then

!Two-Rarefaction Riemann Solver
!p_star
pS=aL+aR-0.5*(gam-1.)*(uR-uL)
pS=pS/(aL/(pL**gamee)+aR/(pR**gamee))
pS=pS**(1./gamee)
!u_star
uS=uL-2.*aL*gamul*(-1.+(pS/pL)**gamee)
!rho_star
rhoSL=rhoL*(pS/pL)**(1./gam)
rhoSR=rhoR*(pS/pR)**(1./gam)

else if(pS>pMin) then

!Two-Shock Riemann Solver
!p_star
p0=max(0._8,pS)
pS=gL(p0)*pL+gR(p0)*pR-(uR-uL)
pS=pS/(gL(p0)+gR(p0))
!u_star
uS=0.5*(uL+uR)+0.5*((pS-pR)*gR(p0)-(pS-pL)*gL(p0))
!rho_star
rhoSL=rhoL*((pS/pL)+gamel)/((pS/pL)*gamel+1.)
rhoSR=rhoR*((pS/pR)+gamel)/((pS/pR)*gamel+1.)
end if

aSL=sqrt(gam*pS/rhoSL)
aSR=sqrt(gam*pS/rhoSR)

end subroutine starRegion

function gL(x) result(pxx)
!Input Variables
real*8,intent(in)::x
!Output Variables
real*8::pxx

  pxx=sqrt((2.*gamil/rhoL)/(x+gamel*pL))
end function gL

function gR(x) result(pxx)
!Input Variables
real*8,intent(in)::x
!Output Variables
real*8::pxx

  pxx=sqrt((2.*gamil/rhoR)/(x+gamel*pR))
end function gR

function limiter(x,y) result(z)
  real*8,intent(in)::x,y
  real*8::z,b

  b=1.

  if(y>0._8)then
    z=max(0.,min(b*x,y),min(x,b*y))
  else 
    z=min(0.,max(b*x,y),max(x,b*y))
  end if

end function limiter


subroutine fluidBound()
integer::k

!outflow boundary on the right and reflecting on the left
do k=1,3
  u2(-1,k)=u2(0,k)
  u2(-2,k)=u2(0,k)
  u2(nx,k)=u2(nx-1,k)
  u2(nx+1,k)=u2(nx-1,k)
end do

u2(-1,2)=-u2(0,2)	
u2(-2,2)=-u2(1,2)
!u2(nx,2)=-u2(nx-1,2)	
!u2(nx+1,2)=-u2(nx-2,2)


go to 1003 
!outflow boundary on the right and left
do k=1,3
  u2(-1,k)=u2(0,k)
  u2(-2,k)=u2(0,k)
  u2(nx,k)=u2(nx-1,k)
  u2(nx+1,k)=u2(nx-1,k)
end do
1003 continue

end subroutine fluidBound

end module RoeSolver_entropyfixed_mod
