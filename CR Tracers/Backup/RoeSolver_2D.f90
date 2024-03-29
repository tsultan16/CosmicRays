!Van Leer MUSCL Scheme (2nd order) for Ieal MHD in 2D
!(Dimensionally split)  
!With Constrained Transport for maintaining Div(B)=0

module RoeSolver2D_mod
use constants_mod
implicit none


real*8::dt,dx,dy !time step and spatial resolution
real*8::gam,gamil,gamul,gamel,gamee,gamuu
real*8::x,y
real*8::rhoL,uL,vL,wL,bxL,byL,bzL,pL,aL !left state
real*8::rhoR,uR,vR,wR,bxR,byR,bzR,pR,aR !right state

!-------------------------------------------------------------
real*8::dens,velx,vely,velz,pres,magx,magy,magz
!-------------------------------------------------------------
real*8::smax,cs,cf,cA
real*8::rhot,ut,vt,wt,pt,byt,bzt,at,ptotL,ptotR,ptot
real*8::lambda(7)
real*8::L1(7),L2(7),L3(7),L4(7),L5(7),L6(7),L7(7)	
real*8::R1(7),R2(7),R3(7),R4(7),R5(7),R6(7),R7(7)
real*8::alpha(7),uright(7),uleft(7)
real*8::alphas,alphaf,betay,betaz,bx,by,bz,bsqr
real*8::BBR,BBL,vvR,vvL,gf1,gf7,gs3,gs5,theta1,theta2,vsq
!-------------------------------------------------------------
real*8::u1(-1:nx+2,-1:ny+2,8),u2(-1:nx+2,-1:ny+2,8),flux(-1:nx+2,8)
real*8::fleft(7),fright(7)
real*8::bxInt(-1:nx+2,-1:ny+2),byInt(-1:nx+2,-1:ny+2)
real*8::fx(-1:nx+2,-1:ny+2),fy(-1:nx+2,-1:ny+2),Ez(-1:nx+2,-1:ny+2)
!-------------------------------------------------------------
character(len=40) :: filename1, filename2, filename3

contains 

!------------------------------------------------------
!initialization
!------------------------------------------------------
subroutine init()

!Local Variables
integer::jj,kk,mode
real*8::vparL,vparR,vperpL,vperpR,bparL,bparR,bperpL,bperpR,kz
real*8::p3,rho3,p1,p2
integer::cell_count(-1:nx+2,-1:ny+2)

open(unit=12,file='Output/horcut_fluid.txt')
!open(unit=13,file='Output/horcut_B.txt')
open(unit=14,file='Output/vercut_fluid.txt')
open(unit=15,file='Output/vercut_B.txt')
open(unit=16,file='Output/diagcut_fluid.txt')
open(unit=17,file='Output/diagcut_B.txt')
open(unit=0+1000,file='Output/t=0.txt')

print*,'First Order Roe Upwind Method...'

!Discountinuity in initial condition at the middle of [xmin,xmax]


dx=(xmax-xmin)/nx
dy=(ymax-ymin)/ny
dt=0.

gam=5._8/3._8
gamil=1./(gam+1.)
gamul=1./(gam-1.)
gamel=(gam-1.)/(gam+1.)
gamee=(gam-1.)/(2.*gam)
gamuu=(gam+1.)/(2.*gam)

!------------------------------------------------------
!Set left and right states 
!------------------------------------------------------
if(initOption==0)then

bxL=1.0
bxR=bxL
byL=0.0
byR=byL

rhoL=1.0
uL=0.
vL=0.
wL=0.
bzL=0.
pL=20.

rhoR=0.125
uR=0.
vR=0.
wR=0.
bzR=0.
pR=5.
end if

if(initOption==1)then

bxL=5./sqrt(4.*pi)
bxR=bxL

rhoL=1.
uL=10.
vL=0.
wL=0.
byL=5./sqrt(4.*pi)
bzL=0.
pL=20.

rhoR=1.
uR=-10.
vR=0.
wR=0.
byR=5./sqrt(4.*pi)
bzR=0.
pR=1.
end if

if(initOption==2)then
bxL=3./sqrt(4.*pi)
bxR=bxL

rhoL=1.0
uL=0.
vL=0.
wL=0.
byL=5./sqrt(4.*pi)
bzL=0.
pL=1.

rhoR=0.1
uR=0.
vR=0.
wR=0.
byR=2./sqrt(4.*pi)
bzR=0.
pR=10.
end if

if(initOption==3)then
bxL=2./sqrt(4.*pi)
bxR=bxL

rhoL=1.08
uL=1.2
vL=0.01
wL=0.5
byL=3.6/sqrt(4.*pi)
bzL=2./sqrt(4.*pi)
pL=0.95

rhoR=1.
uR=0.
vR=0.
wR=0.
byR=4./sqrt(4.*pi)
bzR=2./sqrt(4.*pi)
pR=1.
end if

if(initOption==4)then
bxL=3./sqrt(4.*pi)
bxR=bxL

rhoL=1.
uL=0.
vL=0.
wL=0.
byL=6./sqrt(4.*pi)
bzL=0.
pL=1.

rhoR=0.1
uR=0.
vR=2.
wR=1.
byR=1./sqrt(4.*pi)
bzR=0.
pR=10.
end if


if(initOption==5)then
bxL=0.
bxR=bxL

rhoL=0.1
uL=50.
vL=0.
wL=0.
byL=-1./sqrt(4.*pi)
bzL=-2./sqrt(4.*pi)
pL=0.4

rhoR=0.1
uR=0.
vR=0.
wR=0.
byR=1./sqrt(4.*pi)
bzR=2./sqrt(4.*pi)
pR=0.2
end if

if(initOption==6)then
bxL=0.
bxR=bxL

rhoL=1.
uL=-1.
vL=0.
wL=0.
byL=1.
bzL=0.
pL=1.

rhoR=1.
uR=1.
vR=0.
wR=0.
byR=1.
bzR=0.
pR=1.
end if

if(initOption==7)then
bxL=1.
bxR=bxL

rhoL=1.
uL=0.
vL=0.
wL=0.
byL=1.
bzL=0.
pL=1.

rhoR=0.2
uR=0.
vR=0.
wR=0.
byR=0.
bzR=0.
pR=0.1
end if


if(initOption==8)then
bxL=1.3
bxR=bxL

rhoL=0.4
uL=-0.66991
vL=0.98263
wL=0.
byL=0.0025293
bzL=0.
pL=0.52467

rhoR=1.
uR=0.
vR=0.
wR=0.
byR=1.
bzR=0.
pR=1.
end if

if(initOption==9)then
bxL=0.75
bxR=bxL

rhoL=0.65
uL=0.667
vL=-0.257
wL=0.
byL=0.55
bzL=0.
pL=0.5

rhoR=1.
uR=0.4
vR=-0.94
wR=0.
byR=0.
bzR=0.
pR=0.75
end if

if(initOption==10)then
bxL=0.7
bxR=bxL

rhoL=1.
uL=0.
vL=0.
wL=0.
byL=0.
bzL=0.
pL=1.

rhoR=0.3
uR=0.
vR=0.
wR=1.
byR=1.
bzR=0.
pR=0.2
end if


if(initOption==11)then
bxL=0.75
bxR=bxL

rhoL=1.
uL=0.
vL=0.
wL=0.
byL=1.
bzL=0.
pL=1.

rhoR=0.125
uR=0.
vR=0.
wR=0.
byR=-1.
bzR=0.
pR=0.1
end if


if(initOption==12)then
bxL=1.3
bxR=bxL

rhoL=1.
uL=0.
vL=0.
wL=0.
byL=1.
bzL=0.
pL=1.

rhoR=4.0
uR=0.
vR=0.
wR=0.
byR=-1.
bzR=0.
pR=0.4

end if



if(discontinuityDir==1)then
 do jj=-1,nx+2
 do kk=-1,ny+2
  if(jj<nx/2) then
    u1(jj,kk,1)=rhoL
    u1(jj,kk,2)=rhoL*uL
    u1(jj,kk,3)=rhoL*vL
    u1(jj,kk,4)=rhoL*wL
    u1(jj,kk,5)=bxL
    u1(jj,kk,6)=byL
    u1(jj,kk,7)=bzL  
    u1(jj,kk,8)=0.5*rhoL*(uL*uL+vL*vL+wL*wL) &
           +gamul*pL+0.5*(bxL*bxL+byL*byL+bzL*bzL)

    bxInt(jj,kk)=bxL
    byInt(jj,kk)=byL
  else
    u1(jj,kk,1)=rhoR
    u1(jj,kk,2)=rhoR*uR
    u1(jj,kk,3)=rhoR*vR
    u1(jj,kk,4)=rhoR*wR
    u1(jj,kk,5)=bxR
    u1(jj,kk,6)=byR
    u1(jj,kk,7)=bzR  
    u1(jj,kk,8)=0.5*rhoR*(uR*uR+vR*vR+wR*wR) &
           +gamul*pR+0.5*(bxR*bxR+byR*byR+bzR*bzR) 
    bxInt(jj,kk)=bxR 
    byInt(jj,kk)=byR 
  end if
 end do
 end do

else if(discontinuityDir==2)then
 do jj=-1,nx+2
 do kk=-1,ny+2
  if(kk<ny/2) then
    u1(jj,kk,1)=rhoL
    u1(jj,kk,2)=rhoL*vL
    u1(jj,kk,3)=rhoL*uL
    u1(jj,kk,4)=rhoL*wL
    u1(jj,kk,5)=byL
    u1(jj,kk,6)=bxL
    u1(jj,kk,7)=bzL  
    u1(jj,kk,8)=0.5*rhoL*(uL*uL+vL*vL+wL*wL) &
           +gamul*pL+0.5*(bxL*bxL+byL*byL+bzL*bzL)
    bxInt(jj,kk)=byL
    byInt(jj,kk)=bxL
  else
    u1(jj,kk,1)=rhoR
    u1(jj,kk,2)=rhoR*vR
    u1(jj,kk,3)=rhoR*uR
    u1(jj,kk,4)=rhoR*wR
    u1(jj,kk,5)=byR
    u1(jj,kk,6)=bxR
    u1(jj,kk,7)=bzR  
    u1(jj,kk,8)=0.5*rhoR*(uR*uR+vR*vR+wR*wR) &
           +gamul*pR+0.5*(bxR*bxR+byR*byR+bzR*bzR) 
    bxInt(jj,kk)=byR
    byInt(jj,kk)=bxR
  end if
 end do
 end do
else if(discontinuityDir==3)then

 vparL=uL
 vperpL=vL
 bparL=bxL
 bperpL=byL
 vparR=uR
 vperpR=vR
 bparR=bxR
 bperpR=byR

 uL=(vparL-vperpL)/sqrt(2.)
 vL=(vparL+vperpL)/sqrt(2.)
 bxL=(bparL-bperpL)/sqrt(2.)
 byL=(bparL+bperpL)/sqrt(2.)
 uR=(vparR-vperpR)/sqrt(2.)
 vR=(vparR+vperpR)/sqrt(2.)
 bxR=(bparR-bperpR)/sqrt(2.)
 byR=(bparR+bperpR)/sqrt(2.)


 do jj=-1,nx+2
 do kk=-1,ny+2
  if(kk<=nx/2-jj) then
    u1(jj,kk,1)=rhoL
    u1(jj,kk,2)=rhoL*uL
    u1(jj,kk,3)=rhoL*vL
    u1(jj,kk,4)=rhoL*wL
    u1(jj,kk,5)=bxL
    u1(jj,kk,6)=byL
    u1(jj,kk,7)=bzL  
    u1(jj,kk,8)=0.5*rhoL*(uL*uL+vL*vL+wL*wL) &
           +gamul*pL+0.5*(bxL*bxL+byL*byL+bzL*bzL)
    bxInt(jj,kk)=bxL
    byInt(jj,kk)=byL
  else
    u1(jj,kk,1)=rhoR
    u1(jj,kk,2)=rhoR*uR
    u1(jj,kk,3)=rhoR*vR
    u1(jj,kk,4)=rhoR*wR
    u1(jj,kk,5)=bxR
    u1(jj,kk,6)=byR
    u1(jj,kk,7)=bzR  
    u1(jj,kk,8)=0.5*rhoR*(uR*uR+vR*vR+wR*wR) &
           +gamul*pR+0.5*(bxR*bxR+byR*byR+bzR*bzR) 
    bxInt(jj,kk)=bxR
    byInt(jj,kk)=byR
  end if
 end do
 end do

else if(discontinuityDir==4)then
 do jj=-1,nx+2
 do kk=-1,ny+2
  if((jj-nx/2)*(jj-nx/2)+(kk-ny/2)*(kk-ny/2)<=(nx/4)*(nx/4)) then
    u1(jj,kk,1)=rhoL
    u1(jj,kk,2)=rhoL*uL
    u1(jj,kk,3)=rhoL*vL
    u1(jj,kk,4)=rhoL*wL
    u1(jj,kk,5)=bxL
    u1(jj,kk,6)=byL
    u1(jj,kk,7)=bzL  
    u1(jj,kk,8)=0.5*rhoL*(uL*uL+vL*vL+wL*wL) &
           +gamul*pL+0.5*(bxL*bxL+byL*byL+bzL*bzL)
    bxInt(jj,kk)=bxL
    byInt(jj,kk)=byL
  else
    u1(jj,kk,1)=rhoR
    u1(jj,kk,2)=rhoR*uR
    u1(jj,kk,3)=rhoR*vR
    u1(jj,kk,4)=rhoR*wR
    u1(jj,kk,5)=bxR
    u1(jj,kk,6)=byR
    u1(jj,kk,7)=bzR  
    u1(jj,kk,8)=0.5*rhoR*(uR*uR+vR*vR+wR*wR) &
           +gamul*pR+0.5*(bxR*bxR+byR*byR+bzR*bzR) 

    bxInt(jj,kk)=bxR
    byInt(jj,kk)=byR
  end if
 end do
 end do

end if

if(KH_test==1)then
!Kelvin-Helmholtz Instability 

rhoL=1.0
uL=-1.5!-0.75
vL=0.0
bxL=0.01
byL=0.
bzL=0.
wL=0.
pL=10.0

rhoR=3.0
uR=1.5!0.75
vR=0.0
wR=0.
bxR=bxL
byR=0.
bzR=0.
pR=10.0

mode=4

kz=2.*pi*mode/(xmax-xmin)

do jj=-1,nx+2
 x=xmin+(jj-0.5)*dx
 do kk=-1,ny+2
  y=ymin+(kk-0.5)*dy

  !Perturbation in velocity y-component only
  !Exponential term localizes the perturbation near the interface
  
  !Sinusoidal Perturbation
  if(perturbationType==1)then
    vL=0.05*(ymax-ymin)*sin(kz*x)*exp(-300.*(kk-ny/2)**2.)
  !Random Perturbation
  else if(perturbationType==2)then
    vL=0.05*(ymax-ymin)*rand(0)*exp(-300.*(kk-ny/2)**2.)
  end if
  vR=vL
  if(kk<ny/2) then !Below Interface
    u1(jj,kk,1)=rhoL
    u1(jj,kk,2)=rhoL*uL
    u1(jj,kk,3)=rhoL*vL
    u1(jj,kk,4)=rhoL*wL
    u1(jj,kk,5)=bxL
    u1(jj,kk,6)=byL
    u1(jj,kk,7)=bzL  
    u1(jj,kk,8)=0.5*rhoL*(uL*uL+vL*vL+wL*wL) &
           +gamul*pL+0.5*(bxL*bxL+byL*byL+bzL*bzL)
    bxInt(jj,kk)=bxL
    byInt(jj,kk)=byL
  else
    u1(jj,kk,1)=rhoR
    u1(jj,kk,2)=rhoR*uR
    u1(jj,kk,3)=rhoR*vR
    u1(jj,kk,4)=rhoR*wR
    u1(jj,kk,5)=bxR
    u1(jj,kk,6)=byR
    u1(jj,kk,7)=bzR  
    u1(jj,kk,8)=0.5*rhoR*(uR*uR+vL*vL+wR*wR) &
           +gamul*pR+0.5*(bxR*bxR+byR*byR+bzR*bzR) 
    bxInt(jj,kk)=bxR
    byInt(jj,kk)=byR
  end if
 end do   
end do
end if

if(cloudShock==1)then

!Set up uniform low density static background
!Set up circular region at the center of the grid with high density (aka "the cloud")
!Set up cloud to be in pressure balance

bxL=0.0!1.0
bxR=bxL
byL=0.0
byR=byL

rhoL=100.0
uL=0.
vL=0.
wL=0.
bzL=0.
pL=5.

rhoR=1.0
uR=0.
vR=0.
wR=0.
bzR=0.
pR=5.

p3=80.
rho3=4.

do jj=-1,nx+2
 do kk=-1,ny+2
  if((jj-nx/2)*(jj-nx/2)+(kk-ny/2)*(kk-ny/2)<=(nx/10)*(nx/10)) then
    u1(jj,kk,1)=rhoL
    u1(jj,kk,2)=rhoL*uL
    u1(jj,kk,3)=rhoL*vL
    u1(jj,kk,4)=rhoL*wL
    u1(jj,kk,5)=bxL
    u1(jj,kk,6)=byL
    u1(jj,kk,7)=bzL  
    u1(jj,kk,8)=0.5*rhoL*(uL*uL+vL*vL+wL*wL) &
           +gamul*pL+0.5*(bxL*bxL+byL*byL+bzL*bzL)
    bxInt(jj,kk)=bxL
    byInt(jj,kk)=byL
  else
    u1(jj,kk,1)=rhoR
    u1(jj,kk,2)=rhoR*uR
    u1(jj,kk,3)=rhoR*vR
    u1(jj,kk,4)=rhoR*wR
    u1(jj,kk,5)=bxR
    u1(jj,kk,6)=byR
    u1(jj,kk,7)=bzR  
    u1(jj,kk,8)=0.5*rhoR*(uR*uR+vR*vR+wR*wR) &
           +gamul*pR+0.5*(bxR*bxR+byR*byR+bzR*bzR) 

    bxInt(jj,kk)=bxR
    byInt(jj,kk)=byR
  end if
 end do
end do

!Planar Shock entering from the left
do jj=-1,nx+2
 do kk=-1,ny+2
  if(jj<nx/8) then
    u1(jj,kk,1)=rho3
    u1(jj,kk,2)=rho3*uR
    u1(jj,kk,3)=rho3*vR
    u1(jj,kk,4)=rho3*wR
    u1(jj,kk,5)=bxR
    u1(jj,kk,6)=byR
    u1(jj,kk,7)=bzR  
    u1(jj,kk,8)=0.5*rho3*(uR*uR+vR*vR+wR*wR) &
           +gamul*p3+0.5*(bxR*bxR+byR*byR+bzR*bzR) 

    bxInt(jj,kk)=bxR
    byInt(jj,kk)=byR
  end if
 end do
end do

end if


if(RT_test==1)then
!Rayleigh-Taylor Instability 

!Set gravitational uniform field 
gx=0.
gy=-1.0

rhoL=1.0
uL=0.0!-0.75
vL=0.0
bxL=0.0
byL=0.
bzL=0.
wL=0.
pL=10.0

rhoR=5.0
uR=0.0
vR=0.0
wR=0.
bxR=bxL
byR=0.
bzR=0.
pR=pL

mode=6

kz=2.*pi*mode/(xmax-xmin)

do jj=-1,nx+2
 x=xmin+(jj-0.5)*dx
 do kk=-1,ny+2
  y=ymin+(kk-0.5)*dy

  !Perturbation in velocity y-component only
  !Exponential term localizes the perturbation near the interface
  
  !Sinusoidal Perturbation
  if(perturbationType==1)then
    vL=-0.05*(ymax-ymin)*sin(kz*x)*exp(-2.*pi*abs(kk-0.5*ny)*real(mode)/real(nx))
  !Random Perturbation
  else if(perturbationType==2)then
    vL=-0.05*(ymax-ymin)*rand(0)*exp(-2.*pi*abs(kk-0.5*ny)*mode/nx)
  end if
  vR=vL
  if(kk<0.5*ny) then !Below Interface
    u1(jj,kk,1)=rhoL
    u1(jj,kk,2)=rhoL*uL
    u1(jj,kk,3)=rhoL*vL
    u1(jj,kk,4)=rhoL*wL
    u1(jj,kk,5)=bxL
    u1(jj,kk,6)=byL
    u1(jj,kk,7)=bzL 
    p1=pL!-gy*(rhoR*(0.5*ny-0.5)+rhoL*(0.5*ny-kk+0.5))*dy 
    u1(jj,kk,8)=0.5*rhoL*(uL*uL+vL*vL+wL*wL) &
           +gamul*p1+0.5*(bxL*bxL+byL*byL+bzL*bzL)
    bxInt(jj,kk)=bxL
    byInt(jj,kk)=byL
  else
    u1(jj,kk,1)=rhoR
    u1(jj,kk,2)=rhoR*uR
    u1(jj,kk,3)=rhoR*vR
    u1(jj,kk,4)=rhoR*wR
    u1(jj,kk,5)=bxR
    u1(jj,kk,6)=byR
    u1(jj,kk,7)=bzR 
    p2=pR!+gy*rhoR*(kk-ny-0.5)*dy 
    u1(jj,kk,8)=0.5*rhoR*(uR*uR+vL*vL+wR*wR) &
           +gamul*p2+0.5*(bxR*bxR+byR*byR+bzR*bzR) 
    bxInt(jj,kk)=bxR
    byInt(jj,kk)=byR
  end if
 end do   
end do
end if

if(RT_test==2)then
!Circular Rayleigh-Taylor Instability 

rhoL=1.0
uL=0.0!-0.75
vL=0.0
bxL=0.0
byL=0.
bzL=0.
wL=0.
pL=10.0

rhoR=5.0
uR=0.0
vR=0.0
wR=0.
bxR=bxL
byR=0.
bzR=0.
pR=pL

mode=6

kz=2.*pi*mode/(xmax-xmin)

do jj=-1,nx+2
 x=xmin+(jj-0.5)*dx
 do kk=-1,ny+2
  y=ymin+(kk-0.5)*dy

  !Perturbation in velocity y-component only
  !Exponential term localizes the perturbation near the interface
  
  !Sinusoidal Perturbation
  if(perturbationType==1)then
    vL=-0.05*(ymax-ymin)*sin(kz*x)*exp(-2.*pi*abs(kk-0.5*ny)*real(mode)/real(nx))
  !Random Perturbation
  else if(perturbationType==2)then
    vL=-0.05*(ymax-ymin)*rand(0)*exp(-2.*pi*abs(kk-0.5*ny)*mode/nx)
  end if
   vL=0.
   vR=vL
  if((jj-nx/2)*(jj-nx/2)+(kk-ny/2)*(kk-ny/2)<=(nx/3)*(nx/3)) then !Below Interface
    u1(jj,kk,1)=rhoL
    u1(jj,kk,2)=rhoL*uL
    u1(jj,kk,3)=rhoL*vL
    u1(jj,kk,4)=rhoL*wL
    u1(jj,kk,5)=bxL
    u1(jj,kk,6)=byL
    u1(jj,kk,7)=bzL 
    p1=pL!-gy*(rhoR*(0.5*ny-0.5)+rhoL*(0.5*ny-kk+0.5))*dy 
    u1(jj,kk,8)=0.5*rhoL*(uL*uL+vL*vL+wL*wL) &
           +gamul*p1+0.5*(bxL*bxL+byL*byL+bzL*bzL)
    bxInt(jj,kk)=bxL
    byInt(jj,kk)=byL
  else 
    u1(jj,kk,1)=rhoR
    u1(jj,kk,2)=rhoR*uR
    u1(jj,kk,3)=rhoR*vR
    u1(jj,kk,4)=rhoR*wR
    u1(jj,kk,5)=bxR
    u1(jj,kk,6)=byR
    u1(jj,kk,7)=bzR 
    p2=pR!+gy*rhoR*(kk-ny-0.5)*dy 
    u1(jj,kk,8)=0.5*rhoR*(uR*uR+vL*vL+wR*wR) &
           +gamul*p2+0.5*(bxR*bxR+byR*byR+bzR*bzR) 
    bxInt(jj,kk)=bxR
    byInt(jj,kk)=byR
  end if 
  if((jj-nx/2)*(jj-nx/2)+(kk-ny/2)*(kk-ny/2)<(0.334*nx)*(0.334*nx) &
   .and. (jj-nx/2)*(jj-nx/2)+(kk-ny/2)*(kk-ny/2)>(0.332*nx)*(0.332*nx) ) then
    u1(jj,kk,2)=0.01*rhoL*rand(0)
    u1(jj,kk,3)=0.01*rhoL*rand(0)
  end if
 end do   
end do
end if


u2=u1

!-------------------------------
call bound()
u1=u2
!-------------------------------


print*,'Riemann Problem'
print*,'Left State(rho,vx,vy,vz,Bx,By,Bz,p)=',rhoL,uL,vL,wL,bxL,byL,bzL,pL
print*,'Right State(rho,vx,vy,vz,Bx,By,Bz,p)=',rhoR,uR,vR,wR,bxR,byR,bzR,pR

!call vercut()
call horcut()
!call diagcut()

cell_count=0

call fileOutput(0+1000,cell_count)
!close(unit=0+1000)
print*,'Fluid initialization complete...'


end subroutine init


subroutine computeFluxRoe(n,ind,dir)

!Input Variables
integer::n,ind,dir,jj,kk

!Local Variables
real*8::bxLR,LdotR(7),bxx,byy,bzz,att
real*8::qL(-1:nx+2,8),qR(-1:nx+2,8),slope(-1:nx+2,8)
real*8::ubL(-1:nx+2,8),ubR(-1:nx+2,8),FL(8),FR(8)

!--------------------------------------------------
!Limited Slopes: delta_j
!--------------------------------------------------
if(dir==1)then

do jj=0,n+1
  do kk=1,8
    slope(jj,kk)=limiter(u1(jj,ind,kk)-u1(jj-1,ind,kk),u1(jj+1,ind,kk)-u1(jj,ind,kk))
  end do
end do

else if(dir==2)then
do jj=0,n+1
  do kk=1,8
    slope(jj,kk)=limiter(u1(ind,jj,kk)-u1(ind,jj-1,kk),u1(ind,jj+1,kk)-u1(ind,jj,kk))
  end do
end do

end if

!----------------------------------------------
!Boundary Extrapolated Values: uL_i, uR_i
!----------------------------------------------
do jj=0,n+1
  do kk=1,8
   if(dir==1)then
     qL(jj,kk)=u1(jj,ind,kk)-0.5*slope(jj,kk)
     qR(jj,kk)=u1(jj,ind,kk)+0.5*slope(jj,kk)
   else if(dir==2)then
     qL(jj,kk)=u1(ind,jj,kk)-0.5*slope(jj,kk)
     qR(jj,kk)=u1(ind,jj,kk)+0.5*slope(jj,kk)
   end if
  end do
end do


!-------------------------------------------------------------------
!Evolve Boundary Extrapolated values by a half-step: ubarL_i, ubarR_i
!-------------------------------------------------------------------

if(dir==1)then
do jj=0,n+1 
  rhoL=qL(jj,1)
  uL=qL(jj,2)/qL(jj,1)
  vL=qL(jj,3)/qL(jj,1)
  wL=qL(jj,4)/qL(jj,1) 
  bxL=qL(jj,5)
  byL=qL(jj,6)
  bzL=qL(jj,7)
  BBL=bxL*bxL+byL*byL+bzL*bzL
  vvL=uL*uL+vL*vL+wL*wL
  pL=(gam-1.)*(qL(jj,8)-0.5*rhoL*vvL-0.5*BBL)
  aL=sqrt(gam*pL/rhoL)
  ptotL=pL+0.5*BBL

  FL(1)=rhoL*uL
  FL(2)=rhoL*uL*uL+ptotL-bxL*bxL
  FL(3)=rhoL*uL*vL-bxL*byL
  FL(4)=rhoL*uL*wL-bxL*bzL
  FL(5)=0.
  FL(6)=uL*byL-vL*bxL
  FL(7)=uL*bzL-wL*bxL
  FL(8)=uL*(qL(jj,8)+ptotL) &
      -bxL*(bxL*uL+byL*vL+bzL*wL)

  rhoR=qR(jj,1)
  uR=qR(jj,2)/qR(jj,1)
  vR=qR(jj,3)/qR(jj,1)
  wR=qR(jj,4)/qR(jj,1)
  bxR=qR(jj,5)
  byR=qR(jj,6)
  bzR=qR(jj,7)
  BBR=bxR*bxR+byR*byR+bzR*bzR
  vvR=uR*uR+vR*vR+wR*wR
  pR=(gam-1.)*(qR(jj,8)-0.5*rhoR*vvR-0.5*BBR)
  aR=sqrt(gam*pR/rhoR)
  ptotR=pR+0.5*BBR


  FR(1)=rhoR*uR
  FR(2)=rhoR*uR*uR+ptotR-bxR*bxR
  FR(3)=rhoR*uR*vR-bxR*byR
  FR(4)=rhoR*uR*wR-bxR*bzR
  FR(5)=0.
  FR(6)=uR*byR-vR*bxR
  FR(7)=uR*bzR-wR*bxR
  FR(8)=uR*(qR(jj,8)+ptotR) &
      -bxR*(bxR*uR+byR*vR+bzR*wR)

  do kk=1,8
    ubL(jj,kk)=qL(jj,kk)+0.5*(dt/dx)*(FL(kk)-FR(kk))
    ubR(jj,kk)=qR(jj,kk)+0.5*(dt/dx)*(FL(kk)-FR(kk))
  end do

end do

else if(dir==2)then
do jj=0,n+1 
  rhoL=qL(jj,1)
  uL=qL(jj,2)/qL(jj,1)
  vL=qL(jj,3)/qL(jj,1)
  wL=qL(jj,4)/qL(jj,1) 
  bxL=qL(jj,5)
  byL=qL(jj,6)
  bzL=qL(jj,7)
  BBL=bxL*bxL+byL*byL+bzL*bzL
  vvL=uL*uL+vL*vL+wL*wL
  pL=(gam-1.)*(qL(jj,8)-0.5*rhoL*vvL-0.5*BBL)
  aL=sqrt(gam*pL/rhoL)
  ptotL=pL+0.5*BBL

  FL(1)=rhoL*vL
  FL(2)=rhoL*uL*vL-bxL*byL
  FL(3)=rhoL*vL*vL+ptotL-byL*byL
  FL(4)=rhoL*vL*wL-byL*bzL
  FL(5)=vL*bxL-uL*byL
  FL(6)=0.
  FL(7)=vL*bzL-wL*byL
  FL(8)=vL*(qL(jj,8)+ptotL) &
      -byL*(bxL*uL+byL*vL+bzL*wL)

  rhoR=qR(jj,1)
  uR=qR(jj,2)/qR(jj,1)
  vR=qR(jj,3)/qR(jj,1)
  wR=qR(jj,4)/qR(jj,1)
  bxR=qR(jj,5)
  byR=qR(jj,6)
  bzR=qR(jj,7)
  BBR=bxR*bxR+byR*byR+bzR*bzR
  vvR=uR*uR+vR*vR+wR*wR
  pR=(gam-1.)*(qR(jj,8)-0.5*rhoR*vvR-0.5*BBR)
  aR=sqrt(gam*pR/rhoR)
  ptotR=pR+0.5*BBR


  FR(1)=rhoR*vR
  FR(2)=rhoR*vR*uR-bxR*byR
  FR(3)=rhoR*vR*vR+ptotR-byR*byR
  FR(4)=rhoR*vR*wR-byR*bzR
  FR(5)=vR*bxR-uR*byR
  FR(6)=0.
  FR(7)=vR*bzR-wR*byR
  FR(8)=vR*(qR(jj,8)+ptotR) &
      -byR*(bxR*uR+byR*vR+bzR*wR)

  do kk=1,8
    ubL(jj,kk)=qL(jj,kk)+0.5*(dt/dy)*(FL(kk)-FR(kk))
    ubR(jj,kk)=qR(jj,kk)+0.5*(dt/dy)*(FL(kk)-FR(kk))
  end do
end do

end if

!-------------------------------------------------------------------------

do jj=0,n
!---------------------------------
!w_j+1/2(0)
!---------------------------------

!set left and right states for local Riemann problem at cell interface
!x-sweeps
!--------------------------------
if(dir==1)then
!--------------------------------
rhoL=ubR(jj,1)
uL=ubR(jj,2)/ubR(jj,1)
vL=ubR(jj,3)/ubR(jj,1)
wL=ubR(jj,4)/ubR(jj,1)
bxL=ubR(jj,5)
byL=ubR(jj,6)
bzL=ubR(jj,7)

BBL=bxL*bxL+byL*byL+bzL*bzL
vvL=uL*uL+vL*vL+wL*wL
pL=(gam-1.)*(ubR(jj,8)-0.5*rhoL*vvL-0.5*BBL)
aL=sqrt(gam*pL/rhoL)
ptotL=pL+0.5*BBL

uleft(1)=ubR(jj,1)
uleft(2)=ubR(jj,2)
uleft(3)=ubR(jj,3)
uleft(4)=ubR(jj,4)
uleft(5)=ubR(jj,6)
uleft(6)=ubR(jj,7)
uleft(7)=ubR(jj,8)

rhoR=ubL(jj+1,1)
uR=ubL(jj+1,2)/ubL(jj+1,1)
vR=ubL(jj+1,3)/ubL(jj+1,1)
wR=ubL(jj+1,4)/ubL(jj+1,1)
bxR=ubL(jj+1,5)
byR=ubL(jj+1,6)
bzR=ubL(jj+1,7)

BBR=bxR*bxR+byR*byR+bzR*bzR
vvR=uR*uR+vR*vR+wR*wR

pR=(gam-1.)*(ubL(jj+1,8)-0.5*rhoR*vvR-0.5*BBR)
aR=sqrt(gam*pR/rhoR)
ptotR=pR+0.5*BBR

uright(1)=ubL(jj+1,1)
uright(2)=ubL(jj+1,2)
uright(3)=ubL(jj+1,3)
uright(4)=ubL(jj+1,4)
uright(5)=ubL(jj+1,6)
uright(6)=ubL(jj+1,7)
uright(7)=ubL(jj+1,8)

!--------------------------------
else if(dir==2)then
!--------------------------------
rhoL=ubR(jj,1)
uL=ubR(jj,3)/ubR(jj,1)
vL=ubR(jj,2)/ubR(jj,1)
wL=ubR(jj,4)/ubR(jj,1)
bxL=ubR(jj,6)
byL=ubR(jj,5)
bzL=ubR(jj,7)

BBL=bxL*bxL+byL*byL+bzL*bzL
vvL=uL*uL+vL*vL+wL*wL
pL=(gam-1.)*(ubR(jj,8)-0.5*rhoL*vvL-0.5*BBL)
aL=sqrt(gam*pL/rhoL)
ptotL=pL+0.5*BBL

uleft(1)=ubR(jj,1)
uleft(2)=ubR(jj,3)
uleft(3)=ubR(jj,2)
uleft(4)=ubR(jj,4)
uleft(5)=ubR(jj,5)
uleft(6)=ubR(jj,7)
uleft(7)=ubR(jj,8)

rhoR=ubL(jj+1,1)
uR=ubL(jj+1,3)/ubL(jj+1,1)
vR=ubL(jj+1,2)/ubL(jj+1,1)
wR=ubL(jj+1,4)/ubL(jj+1,1)
bxR=ubL(jj+1,6)
byR=ubL(jj+1,5)
bzR=ubL(jj+1,7)

BBR=bxR*bxR+byR*byR+bzR*bzR
vvR=uR*uR+vR*vR+wR*wR

pR=(gam-1.)*(ubL(jj+1,8)-0.5*rhoR*vvR-0.5*BBR)
aR=sqrt(gam*pR/rhoR)
ptotR=pR+0.5*BBR

uright(1)=ubL(jj+1,1)
uright(2)=ubL(jj+1,3)
uright(3)=ubL(jj+1,2)
uright(4)=ubL(jj+1,4)
uright(5)=ubL(jj+1,5)
uright(6)=ubL(jj+1,7)
uright(7)=ubL(jj+1,8)

end if

fleft(1)=rhoL*uL
fleft(2)=rhoL*uL*uL+ptotL-bxL*bxL
fleft(3)=rhoL*uL*vL-bxL*byL
fleft(4)=rhoL*uL*wL-bxL*bzL
fleft(5)=uL*byL-vL*bxL
fleft(6)=uL*bzL-wL*bxL
fleft(7)=uL*(uleft(7)+ptotL) &
      -bxL*(bxL*uL+byL*vL+bzL*wL)


fright(1)=rhoR*uR
fright(2)=rhoR*uR*uR+ptotR-bxR*bxR
fright(3)=rhoR*uR*vR-bxR*byR
fright(4)=rhoR*uR*wR-bxR*bzR
fright(5)=uR*byR-vR*bxR
fright(6)=uR*bzR-wR*bxR
fright(7)=uR*(uright(7)+ptotR) &
      -bxR*(bxR*uR+byR*vR+bzR*wR)

if(debug==1)then
  if(dir==1)then
   print*,'Cell:',jj,ind
  else if(dir==2)then
   print*,'Cell:',ind,jj
  end if
  print*,'Local RP Left State (rho,u,v,w,bx,by,bz,p) =',rhoL,uL,vL,wL,bxL,byL,bzL,pL
  print*,'Local RP Right State (rho,u,v,w,bx,by,bz,p) =',rhoR,uR,vR,wR,bxR,byR,bzR,pR
end if

!-----------------------------
bxLR=bxL
!-----------------------------

!----------------------------------------------
!Roe-Averaged State Variables (arithmetic mean):
!----------------------------------------------
rhot=0.5_8*(rhoL+rhoR)
ut=0.5_8*(uL+uR)
vt=0.5_8*(vL+vR)
wt=0.5_8*(wL+wR)
byt=0.5_8*(byL+byR)
bzt=0.5_8*(bzL+bzR)
ptot=0.5_8*(ptotL+ptotR) !total pressure: fluid+magnetic field
pt=ptot-0.5_8*(bxLR*bxLR+byt*byt+bzt*bzt)
att=gam*pt/rhot
at=sqrt(att)

if(debug==1)then
  print*,'Roe Average State =',rhot,ut,vt,wt,bxLR,byt,bzt,pt
end if

!-------------------------------------------------
!Alfven and magnetosonic wave speeds: c_s,c_A,c_f
!-------------------------------------------------
bxx=(bxLR*bxLR)/rhot
byy=(byt*byt)/rhot
bzz=(bzt*bzt)/rhot
bsqr=bxx+byy+bzz
vsq=ut*ut+vt*vt+wt*wt

if(debug==1)then
 print*,'bx^2,by^2,bz^2=',bxx,byy,bzz
end if

cA=bxx
cA=sqrt(bxx)
cs=0.5_8*(att+bsqr-sqrt((att+bsqr)*(att+bsqr)-4._8*att*bxx))
!cf=0.5_8*(att+bsqr+sqrt((att+bsqr)*(att+bsqr)-4._8*att*bxx))
cf=att+bsqr-cs
cs=sqrt(cs)
cf=sqrt(cf)


if(debug==1)then
  print*,'Slow,Alfven,Fast,Sound wave speeds=',cs,cA,cf,at
end if

!---------------------------------------------
!Characteristic Eigenvalues: lambda_1..7
!---------------------------------------------
lambda(1)=ut-cf
lambda(2)=ut-cA
lambda(3)=ut-cs
lambda(4)=ut
lambda(5)=ut+cs
lambda(6)=ut+cA
lambda(7)=ut+cf

if(debug==1)then
  print*,'Roe average Eigenvalues=',lambda(1),lambda(2),lambda(3) &
        ,lambda(4),lambda(5),lambda(6),lambda(7)
end if

!-------------------------------------------------
!Roe Jacobian Matrix Eigenvectors: R_1..7, L_1..7
!-------------------------------------------------
if(abs(byt)<tol .and. abs(bzt)<tol)then
 betay=1./sqrt(2.)
 betaz=1./sqrt(2.)
else
 betay=byt/sqrt(byt*byt+bzt*bzt)
 betaz=bzt/sqrt(byt*byt+bzt*bzt)
end if

if(debug==1)then
  print*,'betay,betaz=',betay,betaz
end if

if(abs(byt)<tol .and. abs(bzt)<tol .and. abs(cA-at)<tol)then
 alphas=1.
 alphaf=1.
else
 !!***NOTE**** Inside the square roots in renormalization constants, 
 !! can take abs of numerator to avoid negative values caused by numerical precision errors
 alphas=sqrt(abs(cf*cf-at*at)/(cf*cf-cs*cs))
 alphaf=sqrt(abs(cf*cf-cA*cA)/(cf*cf-cs*cs)) 
end if

if(debug==1)then
  print*,'alphas,alphaf=',alphas,alphaf
end if

gf7=gamul*alphaf*cf*cf+alphaf*cf*ut-alphas*cA*sgn(bxLR)*(betay*vt+betaz*wt) &
    +gamul*(gam-2.)*alphaf*(cf*cf-at*at)
gf1=gamul*alphaf*cf*cf-alphaf*cf*ut+alphas*cA*sgn(bxLR)*(betay*vt+betaz*wt) &
    +gamul*(gam-2.)*alphaf*(cf*cf-at*at)
gs5=gamul*alphas*cs*cs+alphas*cs*ut+alphaf*at*sgn(bxLR)*(betay*vt+betaz*wt) &
    +gamul*(gam-2.)*alphas*(cs*cs-at*at)
gs3=gamul*alphas*cs*cs-alphas*cs*ut-alphaf*at*sgn(bxLR)*(betay*vt+betaz*wt) &
    +gamul*(gam-2.)*alphas*(cs*cs-at*at)

theta1=alphaf*alphaf*at*at*(cf*cf+gamul*(2.-gam)*at*at) &
       +alphas*alphas*cf*cf*(cs*cs+gamul*(2.-gam)*at*at)

theta2=(alphaf*alphaf*cf*at+alphas*alphas*cs*cA)*sgn(bxLR)


R1(1)=alphaf
R1(2)=alphaf*(ut-cf)
R1(3)=alphaf*vt+alphas*betay*cA*sgn(bxLR)
R1(4)=alphaf*wt+alphas*betaz*cA*sgn(bxLR)
R1(5)=(alphas*betay*cf)/sqrt(rhot)
R1(6)=(alphas*betaz*cf)/sqrt(rhot)
R1(7)=alphaf*0.5*vsq+gf1

R2(1)=0.
R2(2)=0.
R2(3)=betaz*sgn(bxLR)
R2(4)=-betay*sgn(bxLR)
R2(5)=betaz/sqrt(rhot)
R2(6)=-betay/sqrt(rhot)
R2(7)=(betaz*vt-betay*wt)*sgn(bxLR)

R3(1)=alphas
R3(2)=alphas*(ut-cs)
R3(3)=alphas*vt-alphaf*betay*at*sgn(bxLR)
R3(4)=alphas*wt-alphaf*betaz*at*sgn(bxLR)
R3(5)=-(alphaf*betay*at*at)/(cf*sqrt(rhot))
R3(6)=-(alphaf*betaz*at*at)/(cf*sqrt(rhot))
R3(7)=alphas*0.5*vsq+gs3

R4(1)=1.
R4(2)=ut
R4(3)=vt
R4(4)=wt
R4(5)=0.
R4(6)=0.
R4(7)=0.5*vsq

R5(1)=alphas
R5(2)=alphas*(ut+cs)
R5(3)=alphas*vt+alphaf*betay*at*sgn(bxLR)
R5(4)=alphas*wt+alphaf*betaz*at*sgn(bxLR)
R5(5)=-(alphaf*betay*at*at)/(cf*sqrt(rhot))
R5(6)=-(alphaf*betaz*at*at)/(cf*sqrt(rhot))
R5(7)=alphas*0.5*vsq+gs5

R6(1)=0.
R6(2)=0.
R6(3)=-betaz*sgn(bxLR)
R6(4)=betay*sgn(bxLR)
R6(5)=betaz/sqrt(rhot)
R6(6)=-betay/sqrt(rhot)
R6(7)=-(betaz*vt-betay*wt)*sgn(bxLR)

R7(1)=alphaf
R7(2)=alphaf*(ut+cf)
R7(3)=alphaf*vt-alphas*betay*cA*sgn(bxLR)
R7(4)=alphaf*wt-alphas*betaz*cA*sgn(bxLR)
R7(5)=(alphas*betay*cf)/sqrt(rhot)
R7(6)=(alphas*betaz*cf)/sqrt(rhot)
R7(7)=alphaf*0.5*vsq+gf7


L1(1)=(1./theta1)*0.25*alphaf*at*at*vsq & 
     +(1./theta2)*( 0.5*alphaf*at*ut*sgn(bxLR) &
     -0.5*alphas*cs*(betay*vt+betaz*wt) )

L1(2)=-(1./theta1)*0.5*alphaf*at*at*ut & 
      -(1./theta2)*0.5*alphaf*at*sgn(bxLR)

L1(3)=-(1./theta1)*0.5*alphaf*at*at*vt & 
      +(1./theta2)*0.5*alphas*betay*cs

L1(4)=-(1./theta1)*0.5*alphaf*at*at*wt & 
      +(1./theta2)*0.5*alphas*betaz*cs

L1(5)=(1./theta1)*0.5*alphas*betay*cf* &
      (cs*cs+gamul*(2.-gam)*at*at)*sqrt(rhot)

L1(6)=(1./theta1)*0.5*alphas*betaz*cf* &
      (cs*cs+gamul*(2.-gam)*at*at)*sqrt(rhot)

L1(7)=(1./theta1)*0.5*alphaf*at*at

L2(1)=-0.5*(betaz*vt-betay*wt)*sgn(bxLR)
L2(2)=0.
L2(3)=0.5*betaz*sgn(bxLR)
L2(4)=-0.5*betay*sgn(bxLR)
L2(5)=0.5*betaz*sqrt(rhot)
L2(6)=-0.5*betay*sqrt(rhot)
L2(7)=0.

L3(1)=(1./theta1)*0.25*alphas*cf*cf*vsq & 
     +(1./theta2)*( 0.5*alphas*cA*ut*sgn(bxLR) &
     +0.5*alphaf*cf*(betay*vt+betaz*wt) )

L3(2)=-(1./theta1)*0.5*alphas*cf*cf*ut & 
      -(1./theta2)*0.5*alphas*cA*sgn(bxLR)

L3(3)=-(1./theta1)*0.5*alphas*cf*cf*vt & 
      -(1./theta2)*0.5*alphaf*betay*cf

L3(4)=-(1./theta1)*0.5*alphas*cf*cf*wt & 
      -(1./theta2)*0.5*alphaf*betaz*cf

L3(5)=-(1./theta1)*0.5*alphaf*betay*cf* &
      (cf*cf+gamul*(2.-gam)*at*at)*sqrt(rhot)

L3(6)=-(1./theta1)*0.5*alphaf*betaz*cf* &
      (cf*cf+gamul*(2.-gam)*at*at)*sqrt(rhot)

L3(7)=(1./theta1)*0.5*alphas*cf*cf


L4(1)=1.-(1./theta1)*0.5*vsq*(alphaf*alphaf*at*at &
      +alphas*alphas*cf*cf )

L4(2)=(1./theta1)*ut*(alphaf*alphaf*at*at &
      +alphas*alphas*cf*cf )

L4(3)=(1./theta1)*vt*(alphaf*alphaf*at*at &
      +alphas*alphas*cf*cf )

L4(4)=(1./theta1)*wt*(alphaf*alphaf*at*at &
      +alphas*alphas*cf*cf )


L4(5)=(1./theta1)*alphaf*alphas*betay*cf*(cf*cf-cs*cs)*sqrt(rhot)

L4(6)=(1./theta1)*alphaf*alphas*betaz*cf*(cf*cf-cs*cs)*sqrt(rhot)

L4(7)=-(1./theta1)*(alphaf*alphaf*at*at+alphas*alphas*cf*cf)


L5(1)=(1./theta1)*0.25*alphas*cf*cf*vsq & 
     -(1./theta2)*( 0.5*alphas*cA*ut*sgn(bxLR) &
     +0.5*alphaf*cf*(betay*vt+betaz*wt) )

L5(2)=-(1./theta1)*0.5*alphas*cf*cf*ut & 
      +(1./theta2)*0.5*alphas*cA*sgn(bxLR)

L5(3)=-(1./theta1)*0.5*alphas*cf*cf*vt & 
      +(1./theta2)*0.5*alphaf*betay*cf

L5(4)=-(1./theta1)*0.5*alphas*cf*cf*wt & 
      +(1./theta2)*0.5*alphaf*betaz*cf

L5(5)=-(1./theta1)*0.5*alphaf*betay*cf* &
      (cf*cf+gamul*(2.-gam)*at*at)*sqrt(rhot)

L5(6)=-(1./theta1)*0.5*alphaf*betaz*cf* &
      (cf*cf+gamul*(2.-gam)*at*at)*sqrt(rhot)

L5(7)=(1./theta1)*0.5*alphas*cf*cf


L6(1)=0.5*(betaz*vt-betay*wt)*sgn(bxLR)
L6(2)=0.
L6(3)=-0.5*betaz*sgn(bxLR)
L6(4)=0.5*betay*sgn(bxLR)
L6(5)=0.5*betaz*sqrt(rhot)
L6(6)=-0.5*betay*sqrt(rhot)
L6(7)=0.


L7(1)=(1./theta1)*0.25*alphaf*at*at*vsq & 
     -(1./theta2)*( 0.5*alphaf*at*ut*sgn(bxLR) &
     -0.5*alphas*cs*(betay*vt+betaz*wt) )

L7(2)=-(1./theta1)*0.5*alphaf*at*at*ut & 
      +(1./theta2)*0.5*alphaf*at*sgn(bxLR)

L7(3)=-(1./theta1)*0.5*alphaf*at*at*vt & 
      -(1./theta2)*0.5*alphas*betay*cs

L7(4)=-(1./theta1)*0.5*alphaf*at*at*wt & 
      -(1./theta2)*0.5*alphas*betaz*cs

L7(5)=(1./theta1)*0.5*alphas*betay*cf* &
      (cs*cs+gamul*(2.-gam)*at*at)*sqrt(rhot)

L7(6)=(1./theta1)*0.5*alphas*betaz*cf* &
      (cs*cs+gamul*(2.-gam)*at*at)*sqrt(rhot)

L7(7)=(1./theta1)*0.5*alphaf*at*at

if((at*at-cA*cA)>tol)then
 if(byt<=tol .and. bzt<tol)then
  R3=-1.*R3
  R5=-1.*R5
  L3=-1.*L3
  L5=-1.*L5
 end if
else 
 if(byt<=tol .and. bzt<tol)then
  R1=-1.*R1
  R7=-1.*R7
  L1=-1.*L1
  L7=-1.*L7
 end if
end if

if(dir==1)then
 !print*,'Cell:',jj,ind
else
!print*,'Cell:',ind,jj
end if
!print*,'||L,R||=',LdotR

!---------------------------------------------
!Roe Averaged Wave Strengths:^alpha_1..7
!---------------------------------------------
alpha=0.

do kk=1,7  
   alpha(1)=alpha(1)+L1(kk)*(uright(kk)-uleft(kk))  
   alpha(2)=alpha(2)+L2(kk)*(uright(kk)-uleft(kk))
   alpha(3)=alpha(3)+L3(kk)*(uright(kk)-uleft(kk))
   alpha(4)=alpha(4)+L4(kk)*(uright(kk)-uleft(kk))
   alpha(5)=alpha(5)+L5(kk)*(uright(kk)-uleft(kk))
   alpha(6)=alpha(6)+L6(kk)*(uright(kk)-uleft(kk))
   alpha(7)=alpha(7)+L7(kk)*(uright(kk)-uleft(kk))
end do


!-------------------------------------------------------------
!Numerical Flux: F_j+1/2 = F(v((x/t)=0))
!------------------------------------------------------------- 
do kk=1,7
  flux(jj,kk)=0.5*(fleft(kk)+fright(kk))-0.5*(     &
         abs(lambda(1))*alpha(1)*R1(kk) &
        +abs(lambda(2))*alpha(2)*R2(kk) &
        +abs(lambda(3))*alpha(3)*R3(kk) &  
        +abs(lambda(4))*alpha(4)*R4(kk) &
        +abs(lambda(5))*alpha(5)*R5(kk) &
        +abs(lambda(6))*alpha(6)*R6(kk) &
        +abs(lambda(7))*alpha(7)*R7(kk) )

end do


!-------------------------------------------------------------
!Cell interface electric field
!-------------------------------------------------------------
if(dir==1)then

fx(jj,ind)=0.5*(byL*uL+byR*uR)-0.5*(     &
         abs(lambda(1))*alpha(1)*R1(5) &
        +abs(lambda(2))*alpha(2)*R2(5) &
        +abs(lambda(3))*alpha(3)*R3(5) &  
        +abs(lambda(4))*alpha(4)*R4(5) &
        +abs(lambda(5))*alpha(5)*R5(5) &
        +abs(lambda(6))*alpha(6)*R6(5) &
        +abs(lambda(7))*alpha(7)*R7(5) )

else if(dir==2)then

fy(ind,jj)=0.5*(byL*uL+byR*uR)-0.5*(     &
         abs(lambda(1))*alpha(1)*R1(5) &
        +abs(lambda(2))*alpha(2)*R2(5) &
        +abs(lambda(3))*alpha(3)*R3(5) &  
        +abs(lambda(4))*alpha(4)*R4(5) &
        +abs(lambda(5))*alpha(5)*R5(5) &
        +abs(lambda(6))*alpha(6)*R6(5) &
        +abs(lambda(7))*alpha(7)*R7(5) )
end if


if(debug==1)then
  print*,'flux=',flux(jj,1),flux(jj,2),flux(jj,3), &
                 flux(jj,4),flux(jj,5),flux(jj,6),flux(jj,7)
end if 


end do

end subroutine computeFluxRoe

subroutine computeFluxHLL(n,ind,dir)

!Input Variables
integer::n,ind,dir,jj,kk

!Local Variables
real*8::bxLR,bxx,byy,bzz,att,SL,SR,lambdaL,lambdaR
real*8::qL(-1:nx+2,8),qR(-1:nx+2,8),slope(-1:nx+2,8)
real*8::ubL(-1:nx+2,8),ubR(-1:nx+2,8),FL(8),FR(8)

!--------------------------------------------------
!Limited Slopes: delta_j
!--------------------------------------------------
if(dir==1)then

do jj=0,n+1
  do kk=1,8
    slope(jj,kk)=limiter(u1(jj,ind,kk)-u1(jj-1,ind,kk),u1(jj+1,ind,kk)-u1(jj,ind,kk))
  end do
end do

else if(dir==2)then
do jj=0,n+1
  do kk=1,8
    slope(jj,kk)=limiter(u1(ind,jj,kk)-u1(ind,jj-1,kk),u1(ind,jj+1,kk)-u1(ind,jj,kk))
  end do
end do

end if

!----------------------------------------------
!Boundary Extrapolated Values: uL_i, uR_i
!----------------------------------------------
do jj=0,n+1
  do kk=1,8
   if(dir==1)then
     qL(jj,kk)=u1(jj,ind,kk)-0.5*slope(jj,kk)
     qR(jj,kk)=u1(jj,ind,kk)+0.5*slope(jj,kk)
   else if(dir==2)then
     qL(jj,kk)=u1(ind,jj,kk)-0.5*slope(jj,kk)
     qR(jj,kk)=u1(ind,jj,kk)+0.5*slope(jj,kk)
   end if
  end do
end do

!-------------------------------------------------------------------
!Evolve Boundary Extrapolated values by a half-step: ubarL_i, ubarR_i
!-------------------------------------------------------------------

if(dir==1)then
do jj=0,n+1 
  rhoL=qL(jj,1)
  uL=qL(jj,2)/qL(jj,1)
  vL=qL(jj,3)/qL(jj,1)
  wL=qL(jj,4)/qL(jj,1) 
  bxL=qL(jj,5)
  byL=qL(jj,6)
  bzL=qL(jj,7)
  BBL=bxL*bxL+byL*byL+bzL*bzL
  vvL=uL*uL+vL*vL+wL*wL
  pL=(gam-1.)*(qL(jj,8)-0.5*rhoL*vvL-0.5*BBL)
  aL=sqrt(gam*pL/rhoL)
  ptotL=pL+0.5*BBL

  FL(1)=rhoL*uL
  FL(2)=rhoL*uL*uL+ptotL-bxL*bxL
  FL(3)=rhoL*uL*vL-bxL*byL
  FL(4)=rhoL*uL*wL-bxL*bzL
  FL(5)=0.
  FL(6)=uL*byL-vL*bxL
  FL(7)=uL*bzL-wL*bxL
  FL(8)=uL*(qL(jj,8)+ptotL) &
      -bxL*(bxL*uL+byL*vL+bzL*wL)

  rhoR=qR(jj,1)
  uR=qR(jj,2)/qR(jj,1)
  vR=qR(jj,3)/qR(jj,1)
  wR=qR(jj,4)/qR(jj,1)
  bxR=qR(jj,5)
  byR=qR(jj,6)
  bzR=qR(jj,7)
  BBR=bxR*bxR+byR*byR+bzR*bzR
  vvR=uR*uR+vR*vR+wR*wR
  pR=(gam-1.)*(qR(jj,8)-0.5*rhoR*vvR-0.5*BBR)
  aR=sqrt(gam*pR/rhoR)
  ptotR=pR+0.5*BBR


  FR(1)=rhoR*uR
  FR(2)=rhoR*uR*uR+ptotR-bxR*bxR
  FR(3)=rhoR*uR*vR-bxR*byR
  FR(4)=rhoR*uR*wR-bxR*bzR
  FR(5)=0.
  FR(6)=uR*byR-vR*bxR
  FR(7)=uR*bzR-wR*bxR
  FR(8)=uR*(qR(jj,8)+ptotR) &
      -bxR*(bxR*uR+byR*vR+bzR*wR)

  do kk=1,8
    ubL(jj,kk)=qL(jj,kk)+0.5*(dt/dx)*(FL(kk)-FR(kk))
    ubR(jj,kk)=qR(jj,kk)+0.5*(dt/dx)*(FL(kk)-FR(kk))
  end do

end do

else if(dir==2)then
do jj=0,n+1 
  rhoL=qL(jj,1)
  uL=qL(jj,2)/qL(jj,1)
  vL=qL(jj,3)/qL(jj,1)
  wL=qL(jj,4)/qL(jj,1) 
  bxL=qL(jj,5)
  byL=qL(jj,6)
  bzL=qL(jj,7)
  BBL=bxL*bxL+byL*byL+bzL*bzL
  vvL=uL*uL+vL*vL+wL*wL
  pL=(gam-1.)*(qL(jj,8)-0.5*rhoL*vvL-0.5*BBL)
  aL=sqrt(gam*pL/rhoL)
  ptotL=pL+0.5*BBL

  FL(1)=rhoL*vL
  FL(2)=rhoL*uL*vL-bxL*byL
  FL(3)=rhoL*vL*vL+ptotL-byL*byL
  FL(4)=rhoL*vL*wL-byL*bzL
  FL(5)=vL*bxL-uL*byL
  FL(6)=0.
  FL(7)=vL*bzL-wL*byL
  FL(8)=vL*(qL(jj,8)+ptotL) &
      -byL*(bxL*uL+byL*vL+bzL*wL)

  rhoR=qR(jj,1)
  uR=qR(jj,2)/qR(jj,1)
  vR=qR(jj,3)/qR(jj,1)
  wR=qR(jj,4)/qR(jj,1)
  bxR=qR(jj,5)
  byR=qR(jj,6)
  bzR=qR(jj,7)
  BBR=bxR*bxR+byR*byR+bzR*bzR
  vvR=uR*uR+vR*vR+wR*wR
  pR=(gam-1.)*(qR(jj,8)-0.5*rhoR*vvR-0.5*BBR)
  aR=sqrt(gam*pR/rhoR)
  ptotR=pR+0.5*BBR


  FR(1)=rhoR*vR
  FR(2)=rhoR*vR*uR-bxR*byR
  FR(3)=rhoR*vR*vR+ptotR-byR*byR
  FR(4)=rhoR*vR*wR-byR*bzR
  FR(5)=vR*bxR-uR*byR
  FR(6)=0.
  FR(7)=vR*bzR-wR*byR
  FR(8)=vR*(qR(jj,8)+ptotR) &
      -byR*(bxR*uR+byR*vR+bzR*wR)

  do kk=1,8
    ubL(jj,kk)=qL(jj,kk)+0.5*(dt/dy)*(FL(kk)-FR(kk))
    ubR(jj,kk)=qR(jj,kk)+0.5*(dt/dy)*(FL(kk)-FR(kk))
  end do
end do

end if

!-------------------------------------------------------------------------
do jj=0,n
!---------------------------------
!w_j+1/2(0)
!---------------------------------

!set left and right states for local Riemann problem at cell interface
!x-sweeps
!--------------------------------
if(dir==1)then
!--------------------------------
rhoL=ubR(jj,1)
uL=ubR(jj,2)/ubR(jj,1)
vL=ubR(jj,3)/ubR(jj,1)
wL=ubR(jj,4)/ubR(jj,1)
bxL=ubR(jj,5)
byL=ubR(jj,6)
bzL=ubR(jj,7)

BBL=bxL*bxL+byL*byL+bzL*bzL
vvL=uL*uL+vL*vL+wL*wL
pL=(gam-1.)*(ubR(jj,8)-0.5*rhoL*vvL-0.5*BBL)
aL=sqrt(gam*pL/rhoL)
ptotL=pL+0.5*BBL

uleft(1)=ubR(jj,1)
uleft(2)=ubR(jj,2)
uleft(3)=ubR(jj,3)
uleft(4)=ubR(jj,4)
uleft(5)=ubR(jj,6)
uleft(6)=ubR(jj,7)
uleft(7)=ubR(jj,8)

rhoR=ubL(jj+1,1)
uR=ubL(jj+1,2)/ubL(jj+1,1)
vR=ubL(jj+1,3)/ubL(jj+1,1)
wR=ubL(jj+1,4)/ubL(jj+1,1)
bxR=ubL(jj+1,5)
byR=ubL(jj+1,6)
bzR=ubL(jj+1,7)

BBR=bxR*bxR+byR*byR+bzR*bzR
vvR=uR*uR+vR*vR+wR*wR

pR=(gam-1.)*(ubL(jj+1,8)-0.5*rhoR*vvR-0.5*BBR)
aR=sqrt(gam*pR/rhoR)
ptotR=pR+0.5*BBR

uright(1)=ubL(jj+1,1)
uright(2)=ubL(jj+1,2)
uright(3)=ubL(jj+1,3)
uright(4)=ubL(jj+1,4)
uright(5)=ubL(jj+1,6)
uright(6)=ubL(jj+1,7)
uright(7)=ubL(jj+1,8)

!--------------------------------
else if(dir==2)then
!--------------------------------
rhoL=ubR(jj,1)
uL=ubR(jj,3)/ubR(jj,1)
vL=ubR(jj,2)/ubR(jj,1)
wL=ubR(jj,4)/ubR(jj,1)
bxL=ubR(jj,6)
byL=ubR(jj,5)
bzL=ubR(jj,7)

BBL=bxL*bxL+byL*byL+bzL*bzL
vvL=uL*uL+vL*vL+wL*wL
pL=(gam-1.)*(ubR(jj,8)-0.5*rhoL*vvL-0.5*BBL)
aL=sqrt(gam*pL/rhoL)
ptotL=pL+0.5*BBL

uleft(1)=ubR(jj,1)
uleft(2)=ubR(jj,3)
uleft(3)=ubR(jj,2)
uleft(4)=ubR(jj,4)
uleft(5)=ubR(jj,5)
uleft(6)=ubR(jj,7)
uleft(7)=ubR(jj,8)

rhoR=ubL(jj+1,1)
uR=ubL(jj+1,3)/ubL(jj+1,1)
vR=ubL(jj+1,2)/ubL(jj+1,1)
wR=ubL(jj+1,4)/ubL(jj+1,1)
bxR=ubL(jj+1,6)
byR=ubL(jj+1,5)
bzR=ubL(jj+1,7)

BBR=bxR*bxR+byR*byR+bzR*bzR
vvR=uR*uR+vR*vR+wR*wR

pR=(gam-1.)*(ubL(jj+1,8)-0.5*rhoR*vvR-0.5*BBR)
aR=sqrt(gam*pR/rhoR)
ptotR=pR+0.5*BBR

uright(1)=ubL(jj+1,1)
uright(2)=ubL(jj+1,3)
uright(3)=ubL(jj+1,2)
uright(4)=ubL(jj+1,4)
uright(5)=ubL(jj+1,5)
uright(6)=ubL(jj+1,7)
uright(7)=ubL(jj+1,8)

end if

fleft(1)=rhoL*uL
fleft(2)=rhoL*uL*uL+ptotL-bxL*bxL
fleft(3)=rhoL*uL*vL-bxL*byL
fleft(4)=rhoL*uL*wL-bxL*bzL
fleft(5)=uL*byL !!!Advective flux for magnetic field normal component/EMF
fleft(6)=uL*bzL-wL*bxL
fleft(7)=uL*(uleft(7)+ptotL) &
      -bxL*(bxL*uL+byL*vL+bzL*wL)


fright(1)=rhoR*uR
fright(2)=rhoR*uR*uR+ptotR-bxR*bxR
fright(3)=rhoR*uR*vR-bxR*byR
fright(4)=rhoR*uR*wR-bxR*bzR
fright(5)=uR*byR !!!Advective flux for magnetic field normal component/EMF
fright(6)=uR*bzR-wR*bxR
fright(7)=uR*(uright(7)+ptotR) &
      -bxR*(bxR*uR+byR*vR+bzR*wR)

if(debug==1)then
  if(dir==1)then
   print*,'Cell:',jj,ind
  else if(dir==2)then
   print*,'Cell:',ind,jj
  end if
  print*,'Local RP Left State (rho,u,v,w,bx,by,bz,p) =',rhoL,uL,vL,wL,bxL,byL,bzL,pL
  print*,'Local RP Right State (rho,u,v,w,bx,by,bz,p) =',rhoR,uR,vR,wR,bxR,byR,bzR,pR
end if

!-----------------------------
bxLR=bxL
!-----------------------------

!----------------------------------------------
!Roe-Averaged State Variables (arithmetic mean):
!----------------------------------------------
rhot=0.5_8*(rhoL+rhoR)
ut=0.5_8*(uL+uR)
vt=0.5_8*(vL+vR)
wt=0.5_8*(wL+wR)
byt=0.5_8*(byL+byR)
bzt=0.5_8*(bzL+bzR)
ptot=0.5_8*(ptotL+ptotR) !total pressure: fluid+magnetic field
pt=ptot-0.5_8*(bxLR*bxLR+byt*byt+bzt*bzt)
att=gam*pt/rhot
at=sqrt(att)

if(debug==1)then
  print*,'Roe Average State =',rhot,ut,vt,wt,bxLR,byt,bzt,pt
end if

!-------------------------------------------------
!Alfven and magnetosonic wave speeds: c_s,c_A,c_f
!-------------------------------------------------
bxx=(bxLR*bxLR)/rhot
byy=(byt*byt)/rhot
bzz=(bzt*bzt)/rhot
bsqr=bxx+byy+bzz
vsq=ut*ut+vt*vt+wt*wt

if(debug==1)then
 print*,'bx^2,by^2,bz^2=',bxx,byy,bzz
end if

cA=bxx
cA=sqrt(bxx)
cs=0.5_8*(att+bsqr-sqrt((att+bsqr)*(att+bsqr)-4._8*att*bxx))
cf=0.5_8*(att+bsqr+sqrt((att+bsqr)*(att+bsqr)-4._8*att*bxx))
!cf=att+bsqr-cs
cs=sqrt(cs)
cf=sqrt(cf)


if(debug==1)then
  print*,'Slow,Alfven,Fast,Sound wave speeds=',cs,cA,cf,at
end if


!----------------------------------------------
!Extremal Wave Speeds
!----------------------------------------------

lambdaL=0.5*(aL*aL+(BBL/rhoL)) &
+0.5*sqrt( (aL*aL+(BBL/rhoL))**2.-4._8*aL*aL*bxL*bxL/rhoL)
lambdaL=uL-sqrt(lambdaL)

lambdaR=0.5*(aR*aR+(BBR/rhoR)) &
+0.5*sqrt( (aR*aR+(BBR/rhoR))**2.-4._8*aR*aR*bxR*bxR/rhoR)
lambdaR=uR+sqrt(lambdaR)

if(HLLOption==0)then
 SL=min(lambdaL,ut-cf)
 SR=max(lambdaR,ut+cf)
else if(HLLOption==1)then
 SL=min(lambdaL,ut-cf,-eps)
 SR=max(lambdaR,ut+cf,eps)
end if

!-------------------------------------------------------------
!Numerical Flux: F_j+1/2 = F(v((x/t)=0))
!------------------------------------------------------------- 
do kk=1,7 
  if(SL>=0.)then
    flux(jj,kk)=fleft(kk)
  else if(SR<=0.)then
    flux(jj,kk)=fright(kk)
  else if(SL<0. .and. SR>0.)then
    flux(jj,kk)=SR*fleft(kk)-SL*fright(kk)+SL*SR*(uright(kk)-uleft(kk))
    flux(jj,kk)=flux(jj,kk)/(SR-SL)
  end if
end do

!-------------------------------------------------------------
!Cell interface electric field
!-------------------------------------------------------------
if(dir==1)then
  fx(jj,ind)=flux(jj,5)
else if(dir==2)then
  fy(ind,jj)=flux(jj,5)
end if



if(debug==1)then
  print*,'flux=',flux(jj,1),flux(jj,2),flux(jj,3), &
                 flux(jj,4),flux(jj,5),flux(jj,6),flux(jj,7)
end if 

end do

end subroutine computeFluxHLL

subroutine computeFluxHLLI(n,ind,dir)

!Input Variables
integer::n,ind,dir,jj,kk,ll,mm

!Local Variables
real*8::bxLR,bxx,byy,bzz,att,SL,SR,lambdaL,lambdaR
!real*8::VLx(7),VRx(7),DUDV(7,7),Rx(7,7),Lx(7,7),Rxu(7,7),Lxu(7,7)
real*8::delta(7)
real*8::qL(-1:nx+2,8),qR(-1:nx+2,8),slope(-1:nx+2,8)
real*8::ubL(-1:nx+2,8),ubR(-1:nx+2,8),FL(8),FR(8)

!--------------------------------------------------
!Limited Slopes: delta_j
!--------------------------------------------------
if(dir==1)then

!$omp parallel do private(kk)
do jj=0,n+1
  do kk=1,8
    slope(jj,kk)=limiter(u1(jj,ind,kk)-u1(jj-1,ind,kk),u1(jj+1,ind,kk)-u1(jj,ind,kk))
  end do
end do
!$omp end parallel do 

else if(dir==2)then
!$omp parallel do private(kk)
do jj=0,n+1
  do kk=1,8
    slope(jj,kk)=limiter(u1(ind,jj,kk)-u1(ind,jj-1,kk),u1(ind,jj+1,kk)-u1(ind,jj,kk))
  end do
end do
!$omp end parallel do 
end if

!----------------------------------------------
!Boundary Extrapolated Values: uL_i, uR_i
!----------------------------------------------
!$omp parallel do private(kk)
do jj=0,n+1
  do kk=1,8
   if(dir==1)then
     qL(jj,kk)=u1(jj,ind,kk)-0.5*slope(jj,kk)
     qR(jj,kk)=u1(jj,ind,kk)+0.5*slope(jj,kk)
   else if(dir==2)then
     qL(jj,kk)=u1(ind,jj,kk)-0.5*slope(jj,kk)
     qR(jj,kk)=u1(ind,jj,kk)+0.5*slope(jj,kk)
   end if
  end do
end do
!$omp end parallel do 

!-------------------------------------------------------------------
!Evolve Boundary Extrapolated values by a half-step: ubarL_i, ubarR_i
!-------------------------------------------------------------------

if(dir==1)then
!$omp parallel do private(kk,rhoL,uL,vL,wL,bxL,byL,bzL,BBL,vvL,pL,aL,ptotL, &
!$omp& FL,rhoR,uR,vR,wR,bxR,byR,bzR,BBR,vvR,pR,aR,ptotR,FR)
do jj=0,n+1 
  rhoL=qL(jj,1)
  uL=qL(jj,2)/qL(jj,1)
  vL=qL(jj,3)/qL(jj,1)
  wL=qL(jj,4)/qL(jj,1) 
  bxL=qL(jj,5)
  byL=qL(jj,6)
  bzL=qL(jj,7)
  BBL=bxL*bxL+byL*byL+bzL*bzL
  vvL=uL*uL+vL*vL+wL*wL
  pL=(gam-1.)*(qL(jj,8)-0.5*rhoL*vvL-0.5*BBL)
  aL=sqrt(gam*pL/rhoL)
  ptotL=pL+0.5*BBL

  FL(1)=rhoL*uL
  FL(2)=rhoL*uL*uL+ptotL-bxL*bxL
  FL(3)=rhoL*uL*vL-bxL*byL
  FL(4)=rhoL*uL*wL-bxL*bzL
  FL(5)=0.
  FL(6)=uL*byL-vL*bxL
  FL(7)=uL*bzL-wL*bxL
  FL(8)=uL*(qL(jj,8)+ptotL) &
      -bxL*(bxL*uL+byL*vL+bzL*wL)

  rhoR=qR(jj,1)
  uR=qR(jj,2)/qR(jj,1)
  vR=qR(jj,3)/qR(jj,1)
  wR=qR(jj,4)/qR(jj,1)
  bxR=qR(jj,5)
  byR=qR(jj,6)
  bzR=qR(jj,7)
  BBR=bxR*bxR+byR*byR+bzR*bzR
  vvR=uR*uR+vR*vR+wR*wR
  pR=(gam-1.)*(qR(jj,8)-0.5*rhoR*vvR-0.5*BBR)
  aR=sqrt(gam*pR/rhoR)
  ptotR=pR+0.5*BBR


  FR(1)=rhoR*uR
  FR(2)=rhoR*uR*uR+ptotR-bxR*bxR
  FR(3)=rhoR*uR*vR-bxR*byR
  FR(4)=rhoR*uR*wR-bxR*bzR
  FR(5)=0.
  FR(6)=uR*byR-vR*bxR
  FR(7)=uR*bzR-wR*bxR
  FR(8)=uR*(qR(jj,8)+ptotR) &
      -bxR*(bxR*uR+byR*vR+bzR*wR)

  do kk=1,8
    ubL(jj,kk)=qL(jj,kk)+0.5*(dt/dx)*(FL(kk)-FR(kk))
    ubR(jj,kk)=qR(jj,kk)+0.5*(dt/dx)*(FL(kk)-FR(kk))
  end do

end do
!$omp end parallel do

else if(dir==2)then
!$omp parallel do private(kk,rhoL,uL,vL,wL,bxL,byL,bzL,BBL,vvL,pL,aL, &
!$omp& ptotL,FL,rhoR,uR,vR,wR,bxR,byR,bzR,BBR,vvR,pR,aR,ptotR,FR)
do jj=0,n+1 
  rhoL=qL(jj,1)
  uL=qL(jj,2)/qL(jj,1)
  vL=qL(jj,3)/qL(jj,1)
  wL=qL(jj,4)/qL(jj,1) 
  bxL=qL(jj,5)
  byL=qL(jj,6)
  bzL=qL(jj,7)
  BBL=bxL*bxL+byL*byL+bzL*bzL
  vvL=uL*uL+vL*vL+wL*wL
  pL=(gam-1.)*(qL(jj,8)-0.5*rhoL*vvL-0.5*BBL)
  aL=sqrt(gam*pL/rhoL)
  ptotL=pL+0.5*BBL

  FL(1)=rhoL*vL
  FL(2)=rhoL*uL*vL-bxL*byL
  FL(3)=rhoL*vL*vL+ptotL-byL*byL
  FL(4)=rhoL*vL*wL-byL*bzL
  FL(5)=vL*bxL-uL*byL
  FL(6)=0.
  FL(7)=vL*bzL-wL*byL
  FL(8)=vL*(qL(jj,8)+ptotL) &
      -byL*(bxL*uL+byL*vL+bzL*wL)

  rhoR=qR(jj,1)
  uR=qR(jj,2)/qR(jj,1)
  vR=qR(jj,3)/qR(jj,1)
  wR=qR(jj,4)/qR(jj,1)
  bxR=qR(jj,5)
  byR=qR(jj,6)
  bzR=qR(jj,7)
  BBR=bxR*bxR+byR*byR+bzR*bzR
  vvR=uR*uR+vR*vR+wR*wR
  pR=(gam-1.)*(qR(jj,8)-0.5*rhoR*vvR-0.5*BBR)
  aR=sqrt(gam*pR/rhoR)
  ptotR=pR+0.5*BBR


  FR(1)=rhoR*vR
  FR(2)=rhoR*vR*uR-bxR*byR
  FR(3)=rhoR*vR*vR+ptotR-byR*byR
  FR(4)=rhoR*vR*wR-byR*bzR
  FR(5)=vR*bxR-uR*byR
  FR(6)=0.
  FR(7)=vR*bzR-wR*byR
  FR(8)=vR*(qR(jj,8)+ptotR) &
      -byR*(bxR*uR+byR*vR+bzR*wR)

  do kk=1,8
    ubL(jj,kk)=qL(jj,kk)+0.5*(dt/dy)*(FL(kk)-FR(kk))
    ubR(jj,kk)=qR(jj,kk)+0.5*(dt/dy)*(FL(kk)-FR(kk))
  end do

end do
!$omp end parallel do

end if

!-------------------------------------------------------------------------
!$omp parallel do shared(ubL,ubR,flux,fx,fy) &
!$omp& private(kk,rhoL,uL,vL,wL,bxL,byL,bzL,BBL,vvL,pL, &
!$omp& aL,ptotL,fleft,bxLR,rhoR,uR,vR,wR,bxR,byR,bzR,BBR,vvR,pR, & 
!$omp& aR,ptotR,fright,uleft,uright,rhot,ut,vt,wt,pt,byt,bzt,at, &
!$omp& att,ptot,bxx,byy,bzz,bsqr,vsq,cA,cs,cf,betay, &
!$omp& betaz,alphas,alphaf,gf1,gf7,gs3,gs5,theta1,theta2,lambda, &
!$omp& R1,R2,R3,R4,R5,R6,R7,L1,L2,L3,L4,L5,L6,L7,lambdaL,lambdaR, &
!$omp& SL,SR,delta,alpha)
do jj=0,n
!---------------------------------
!w_j+1/2(0)
!---------------------------------

!set left and right states for local Riemann problem at cell interface
!x-sweeps
!--------------------------------
if(dir==1)then
!--------------------------------
rhoL=ubR(jj,1)
uL=ubR(jj,2)/ubR(jj,1)
vL=ubR(jj,3)/ubR(jj,1)
wL=ubR(jj,4)/ubR(jj,1)
bxL=ubR(jj,5)
byL=ubR(jj,6)
bzL=ubR(jj,7)

BBL=bxL*bxL+byL*byL+bzL*bzL
vvL=uL*uL+vL*vL+wL*wL
pL=(gam-1.)*(ubR(jj,8)-0.5*rhoL*vvL-0.5*BBL)
aL=sqrt(gam*pL/rhoL)
ptotL=pL+0.5*BBL

uleft(1)=ubR(jj,1)
uleft(2)=ubR(jj,2)
uleft(3)=ubR(jj,3)
uleft(4)=ubR(jj,4)
uleft(5)=ubR(jj,6)
uleft(6)=ubR(jj,7)
uleft(7)=ubR(jj,8)

rhoR=ubL(jj+1,1)
uR=ubL(jj+1,2)/ubL(jj+1,1)
vR=ubL(jj+1,3)/ubL(jj+1,1)
wR=ubL(jj+1,4)/ubL(jj+1,1)
bxR=ubL(jj+1,5)
byR=ubL(jj+1,6)
bzR=ubL(jj+1,7)

BBR=bxR*bxR+byR*byR+bzR*bzR
vvR=uR*uR+vR*vR+wR*wR

pR=(gam-1.)*(ubL(jj+1,8)-0.5*rhoR*vvR-0.5*BBR)
aR=sqrt(gam*pR/rhoR)
ptotR=pR+0.5*BBR

uright(1)=ubL(jj+1,1)
uright(2)=ubL(jj+1,2)
uright(3)=ubL(jj+1,3)
uright(4)=ubL(jj+1,4)
uright(5)=ubL(jj+1,6)
uright(6)=ubL(jj+1,7)
uright(7)=ubL(jj+1,8)

!--------------------------------
else if(dir==2)then
!--------------------------------
rhoL=ubR(jj,1)
uL=ubR(jj,3)/ubR(jj,1)
vL=ubR(jj,2)/ubR(jj,1)
wL=ubR(jj,4)/ubR(jj,1)
bxL=ubR(jj,6)
byL=ubR(jj,5)
bzL=ubR(jj,7)

BBL=bxL*bxL+byL*byL+bzL*bzL
vvL=uL*uL+vL*vL+wL*wL
pL=(gam-1.)*(ubR(jj,8)-0.5*rhoL*vvL-0.5*BBL)
aL=sqrt(gam*pL/rhoL)
ptotL=pL+0.5*BBL

uleft(1)=ubR(jj,1)
uleft(2)=ubR(jj,3)
uleft(3)=ubR(jj,2)
uleft(4)=ubR(jj,4)
uleft(5)=ubR(jj,5)
uleft(6)=ubR(jj,7)
uleft(7)=ubR(jj,8)

rhoR=ubL(jj+1,1)
uR=ubL(jj+1,3)/ubL(jj+1,1)
vR=ubL(jj+1,2)/ubL(jj+1,1)
wR=ubL(jj+1,4)/ubL(jj+1,1)
bxR=ubL(jj+1,6)
byR=ubL(jj+1,5)
bzR=ubL(jj+1,7)

BBR=bxR*bxR+byR*byR+bzR*bzR
vvR=uR*uR+vR*vR+wR*wR

pR=(gam-1.)*(ubL(jj+1,8)-0.5*rhoR*vvR-0.5*BBR)
aR=sqrt(gam*pR/rhoR)
ptotR=pR+0.5*BBR

uright(1)=ubL(jj+1,1)
uright(2)=ubL(jj+1,3)
uright(3)=ubL(jj+1,2)
uright(4)=ubL(jj+1,4)
uright(5)=ubL(jj+1,5)
uright(6)=ubL(jj+1,7)
uright(7)=ubL(jj+1,8)

end if

fleft(1)=rhoL*uL
fleft(2)=rhoL*uL*uL+ptotL-bxL*bxL
fleft(3)=rhoL*uL*vL-bxL*byL
fleft(4)=rhoL*uL*wL-bxL*bzL
fleft(5)=uL*byL !!!Advective flux for magnetic field normal component/EMF
fleft(6)=uL*bzL-wL*bxL
fleft(7)=uL*(uleft(7)+ptotL) &
      -bxL*(bxL*uL+byL*vL+bzL*wL)


fright(1)=rhoR*uR
fright(2)=rhoR*uR*uR+ptotR-bxR*bxR
fright(3)=rhoR*uR*vR-bxR*byR
fright(4)=rhoR*uR*wR-bxR*bzR
fright(5)=uR*byR !!!Advective flux for magnetic field normal component/EMF
fright(6)=uR*bzR-wR*bxR
fright(7)=uR*(uright(7)+ptotR) &
      -bxR*(bxR*uR+byR*vR+bzR*wR)


if(debug==1)then
  if(dir==1)then
   print*,'Cell:',jj,ind
  else if(dir==2)then
   print*,'Cell:',ind,jj
  end if
  print*,'Local RP Left State (rho,u,v,w,bx,by,bz,p) =',rhoL,uL,vL,wL,bxL,byL,bzL,pL
  print*,'Local RP Right State (rho,u,v,w,bx,by,bz,p) =',rhoR,uR,vR,wR,bxR,byR,bzR,pR
end if

!-----------------------------
bxLR=bxL
!-----------------------------

!----------------------------------------------
!Roe-Averaged State Variables (arithmetic mean):
!----------------------------------------------
rhot=0.5_8*(rhoL+rhoR)
ut=0.5_8*(uL+uR)
vt=0.5_8*(vL+vR)
wt=0.5_8*(wL+wR)
byt=0.5_8*(byL+byR)
bzt=0.5_8*(bzL+bzR)
ptot=0.5_8*(ptotL+ptotR) !total pressure: fluid+magnetic field
pt=ptot-0.5_8*(bxLR*bxLR+byt*byt+bzt*bzt)
att=gam*pt/rhot
at=sqrt(att)

if(debug==1)then
  print*,'Roe Average State =',rhot,ut,vt,wt,bxLR,byt,bzt,pt
end if

!-------------------------------------------------
!Alfven and magnetosonic wave speeds: c_s,c_A,c_f
!-------------------------------------------------
bxx=(bxLR*bxLR)/rhot
byy=(byt*byt)/rhot
bzz=(bzt*bzt)/rhot
bsqr=bxx+byy+bzz
vsq=ut*ut+vt*vt+wt*wt

if(debug==1)then
 print*,'bx^2,by^2,bz^2=',bxx,byy,bzz
end if

cA=bxx
cA=sqrt(bxx)
cs=0.5_8*(att+bsqr-sqrt((att+bsqr)*(att+bsqr)-4._8*att*bxx))
cf=0.5_8*(att+bsqr+sqrt((att+bsqr)*(att+bsqr)-4._8*att*bxx))
!cf=att+bsqr-cs
cs=sqrt(cs)
cf=sqrt(cf)


if(debug==1)then
  print*,'Slow,Alfven,Fast,Sound wave speeds=',cs,cA,cf,at
end if


!---------------------------------------------
!Roe Averaged Characteristic Eigenvalues: lambda_1..7
!---------------------------------------------
lambda(1)=ut-cf
lambda(2)=ut-cA
lambda(3)=ut-cs
lambda(4)=ut
lambda(5)=ut+cs
lambda(6)=ut+cA
lambda(7)=ut+cf


!-------------------------------------------------
!Roe Jacobian Matrix Eigenvectors: R_1..7, L_1..7
!-------------------------------------------------
if(abs(byt)<tol .and. abs(bzt)<tol)then
 betay=1./sqrt(2.)
 betaz=1./sqrt(2.)
else
 betay=byt/sqrt(byt*byt+bzt*bzt)
 betaz=bzt/sqrt(byt*byt+bzt*bzt)
end if

if(debug==1)then
  print*,'betay,betaz=',betay,betaz
end if

if(abs(byt)<tol .and. abs(bzt)<tol .and. abs(cA-at)<tol)then
 alphas=1.
 alphaf=1.
else
 !!***NOTE**** Inside the square roots in renormalization constants, 
 !! can take abs of numerator to avoid negative values caused by numerical precision errors
 alphas=sqrt(abs(cf*cf-at*at)/(cf*cf-cs*cs))
 alphaf=sqrt(abs(cf*cf-cA*cA)/(cf*cf-cs*cs)) 
end if

if(debug==1)then
  print*,'alphas,alphaf=',alphas,alphaf
end if

gf7=gamul*alphaf*cf*cf+alphaf*cf*ut-alphas*cA*sgn(bxLR)*(betay*vt+betaz*wt) &
    +gamul*(gam-2.)*alphaf*(cf*cf-at*at)
gf1=gamul*alphaf*cf*cf-alphaf*cf*ut+alphas*cA*sgn(bxLR)*(betay*vt+betaz*wt) &
    +gamul*(gam-2.)*alphaf*(cf*cf-at*at)
gs5=gamul*alphas*cs*cs+alphas*cs*ut+alphaf*at*sgn(bxLR)*(betay*vt+betaz*wt) &
    +gamul*(gam-2.)*alphas*(cs*cs-at*at)
gs3=gamul*alphas*cs*cs-alphas*cs*ut-alphaf*at*sgn(bxLR)*(betay*vt+betaz*wt) &
    +gamul*(gam-2.)*alphas*(cs*cs-at*at)

theta1=alphaf*alphaf*at*at*(cf*cf+gamul*(2.-gam)*at*at) &
       +alphas*alphas*cf*cf*(cs*cs+gamul*(2.-gam)*at*at)

theta2=(alphaf*alphaf*cf*at+alphas*alphas*cs*cA)*sgn(bxLR)


R1(1)=alphaf
R1(2)=alphaf*(ut-cf)
R1(3)=alphaf*vt+alphas*betay*cA*sgn(bxLR)
R1(4)=alphaf*wt+alphas*betaz*cA*sgn(bxLR)
R1(5)=(alphas*betay*cf)/sqrt(rhot)
R1(6)=(alphas*betaz*cf)/sqrt(rhot)
R1(7)=alphaf*0.5*vsq+gf1

R2(1)=0.
R2(2)=0.
R2(3)=betaz*sgn(bxLR)
R2(4)=-betay*sgn(bxLR)
R2(5)=betaz/sqrt(rhot)
R2(6)=-betay/sqrt(rhot)
R2(7)=(betaz*vt-betay*wt)*sgn(bxLR)

R3(1)=alphas
R3(2)=alphas*(ut-cs)
R3(3)=alphas*vt-alphaf*betay*at*sgn(bxLR)
R3(4)=alphas*wt-alphaf*betaz*at*sgn(bxLR)
R3(5)=-(alphaf*betay*at*at)/(cf*sqrt(rhot))
R3(6)=-(alphaf*betaz*at*at)/(cf*sqrt(rhot))
R3(7)=alphas*0.5*vsq+gs3

R4(1)=1.
R4(2)=ut
R4(3)=vt
R4(4)=wt
R4(5)=0.
R4(6)=0.
R4(7)=0.5*vsq

R5(1)=alphas
R5(2)=alphas*(ut+cs)
R5(3)=alphas*vt+alphaf*betay*at*sgn(bxLR)
R5(4)=alphas*wt+alphaf*betaz*at*sgn(bxLR)
R5(5)=-(alphaf*betay*at*at)/(cf*sqrt(rhot))
R5(6)=-(alphaf*betaz*at*at)/(cf*sqrt(rhot))
R5(7)=alphas*0.5*vsq+gs5

R6(1)=0.
R6(2)=0.
R6(3)=-betaz*sgn(bxLR)
R6(4)=betay*sgn(bxLR)
R6(5)=betaz/sqrt(rhot)
R6(6)=-betay/sqrt(rhot)
R6(7)=-(betaz*vt-betay*wt)*sgn(bxLR)

R7(1)=alphaf
R7(2)=alphaf*(ut+cf)
R7(3)=alphaf*vt-alphas*betay*cA*sgn(bxLR)
R7(4)=alphaf*wt-alphas*betaz*cA*sgn(bxLR)
R7(5)=(alphas*betay*cf)/sqrt(rhot)
R7(6)=(alphas*betaz*cf)/sqrt(rhot)
R7(7)=alphaf*0.5*vsq+gf7


L1(1)=(1./theta1)*0.25*alphaf*at*at*vsq & 
     +(1./theta2)*( 0.5*alphaf*at*ut*sgn(bxLR) &
     -0.5*alphas*cs*(betay*vt+betaz*wt) )

L1(2)=-(1./theta1)*0.5*alphaf*at*at*ut & 
      -(1./theta2)*0.5*alphaf*at*sgn(bxLR)

L1(3)=-(1./theta1)*0.5*alphaf*at*at*vt & 
      +(1./theta2)*0.5*alphas*betay*cs

L1(4)=-(1./theta1)*0.5*alphaf*at*at*wt & 
      +(1./theta2)*0.5*alphas*betaz*cs

L1(5)=(1./theta1)*0.5*alphas*betay*cf* &
      (cs*cs+gamul*(2.-gam)*at*at)*sqrt(rhot)

L1(6)=(1./theta1)*0.5*alphas*betaz*cf* &
      (cs*cs+gamul*(2.-gam)*at*at)*sqrt(rhot)

L1(7)=(1./theta1)*0.5*alphaf*at*at

L2(1)=-0.5*(betaz*vt-betay*wt)*sgn(bxLR)
L2(2)=0.
L2(3)=0.5*betaz*sgn(bxLR)
L2(4)=-0.5*betay*sgn(bxLR)
L2(5)=0.5*betaz*sqrt(rhot)
L2(6)=-0.5*betay*sqrt(rhot)
L2(7)=0.

L3(1)=(1./theta1)*0.25*alphas*cf*cf*vsq & 
     +(1./theta2)*( 0.5*alphas*cA*ut*sgn(bxLR) &
     +0.5*alphaf*cf*(betay*vt+betaz*wt) )

L3(2)=-(1./theta1)*0.5*alphas*cf*cf*ut & 
      -(1./theta2)*0.5*alphas*cA*sgn(bxLR)

L3(3)=-(1./theta1)*0.5*alphas*cf*cf*vt & 
      -(1./theta2)*0.5*alphaf*betay*cf

L3(4)=-(1./theta1)*0.5*alphas*cf*cf*wt & 
      -(1./theta2)*0.5*alphaf*betaz*cf

L3(5)=-(1./theta1)*0.5*alphaf*betay*cf* &
      (cf*cf+gamul*(2.-gam)*at*at)*sqrt(rhot)

L3(6)=-(1./theta1)*0.5*alphaf*betaz*cf* &
      (cf*cf+gamul*(2.-gam)*at*at)*sqrt(rhot)

L3(7)=(1./theta1)*0.5*alphas*cf*cf


L4(1)=1.-(1./theta1)*0.5*vsq*(alphaf*alphaf*at*at &
      +alphas*alphas*cf*cf )

L4(2)=(1./theta1)*ut*(alphaf*alphaf*at*at &
      +alphas*alphas*cf*cf )

L4(3)=(1./theta1)*vt*(alphaf*alphaf*at*at &
      +alphas*alphas*cf*cf )

L4(4)=(1./theta1)*wt*(alphaf*alphaf*at*at &
      +alphas*alphas*cf*cf )


L4(5)=(1./theta1)*alphaf*alphas*betay*cf*(cf*cf-cs*cs)*sqrt(rhot)

L4(6)=(1./theta1)*alphaf*alphas*betaz*cf*(cf*cf-cs*cs)*sqrt(rhot)

L4(7)=-(1./theta1)*(alphaf*alphaf*at*at+alphas*alphas*cf*cf)


L5(1)=(1./theta1)*0.25*alphas*cf*cf*vsq & 
     -(1./theta2)*( 0.5*alphas*cA*ut*sgn(bxLR) &
     +0.5*alphaf*cf*(betay*vt+betaz*wt) )

L5(2)=-(1./theta1)*0.5*alphas*cf*cf*ut & 
      +(1./theta2)*0.5*alphas*cA*sgn(bxLR)

L5(3)=-(1./theta1)*0.5*alphas*cf*cf*vt & 
      +(1./theta2)*0.5*alphaf*betay*cf

L5(4)=-(1./theta1)*0.5*alphas*cf*cf*wt & 
      +(1./theta2)*0.5*alphaf*betaz*cf

L5(5)=-(1./theta1)*0.5*alphaf*betay*cf* &
      (cf*cf+gamul*(2.-gam)*at*at)*sqrt(rhot)

L5(6)=-(1./theta1)*0.5*alphaf*betaz*cf* &
      (cf*cf+gamul*(2.-gam)*at*at)*sqrt(rhot)

L5(7)=(1./theta1)*0.5*alphas*cf*cf


L6(1)=0.5*(betaz*vt-betay*wt)*sgn(bxLR)
L6(2)=0.
L6(3)=-0.5*betaz*sgn(bxLR)
L6(4)=0.5*betay*sgn(bxLR)
L6(5)=0.5*betaz*sqrt(rhot)
L6(6)=-0.5*betay*sqrt(rhot)
L6(7)=0.


L7(1)=(1./theta1)*0.25*alphaf*at*at*vsq & 
     -(1./theta2)*( 0.5*alphaf*at*ut*sgn(bxLR) &
     -0.5*alphas*cs*(betay*vt+betaz*wt) )

L7(2)=-(1./theta1)*0.5*alphaf*at*at*ut & 
      +(1./theta2)*0.5*alphaf*at*sgn(bxLR)

L7(3)=-(1./theta1)*0.5*alphaf*at*at*vt & 
      -(1./theta2)*0.5*alphas*betay*cs

L7(4)=-(1./theta1)*0.5*alphaf*at*at*wt & 
      -(1./theta2)*0.5*alphas*betaz*cs

L7(5)=(1./theta1)*0.5*alphas*betay*cf* &
      (cs*cs+gamul*(2.-gam)*at*at)*sqrt(rhot)

L7(6)=(1./theta1)*0.5*alphas*betaz*cf* &
      (cs*cs+gamul*(2.-gam)*at*at)*sqrt(rhot)

L7(7)=(1./theta1)*0.5*alphaf*at*at

if((at*at-cA*cA)>tol)then
 if(byt<=tol .and. bzt<tol)then
  R3=-1.*R3
  R5=-1.*R5
  L3=-1.*L3
  L5=-1.*L5
 end if
else 
 if(byt<=tol .and. bzt<tol)then
  R1=-1.*R1
  R7=-1.*R7
  L1=-1.*L1
  L7=-1.*L7
 end if
end if

!---------------------------------------------
!Roe Averaged Wave Strengths:^alpha_1..7
!---------------------------------------------
alpha=0.

do kk=1,7  
   alpha(1)=alpha(1)+L1(kk)*(uright(kk)-uleft(kk))  
   alpha(2)=alpha(2)+L2(kk)*(uright(kk)-uleft(kk))
   alpha(3)=alpha(3)+L3(kk)*(uright(kk)-uleft(kk))
   alpha(4)=alpha(4)+L4(kk)*(uright(kk)-uleft(kk))
   alpha(5)=alpha(5)+L5(kk)*(uright(kk)-uleft(kk))
   alpha(6)=alpha(6)+L6(kk)*(uright(kk)-uleft(kk))
   alpha(7)=alpha(7)+L7(kk)*(uright(kk)-uleft(kk))
end do

!----------------------------------------------
!Extremal Wave Speeds
!----------------------------------------------

lambdaL=0.5*(aL*aL+(BBL/rhoL)) &
+0.5*sqrt( (aL*aL+(BBL/rhoL))**2.-4._8*aL*aL*bxL*bxL/rhoL)
lambdaL=uL-sqrt(lambdaL)

lambdaR=0.5*(aR*aR+(BBR/rhoR)) &
+0.5*sqrt( (aR*aR+(BBR/rhoR))**2.-4._8*aR*aR*bxR*bxR/rhoR)
lambdaR=uR+sqrt(lambdaR)

if(HLLOption==0)then
 SL=min(lambdaL,ut-cf)
 SR=max(lambdaR,ut+cf)
else if(HLLOption==1)then
 SL=min(lambdaL,ut-cf,-eps)
 SR=max(lambdaR,ut+cf,eps)
end if

!-------------------------------------------------------------
!Numerical Flux: F_j+1/2 = F(v((x/t)=0))
!------------------------------------------------------------- 
do kk=1,7
  delta(kk)=1.-(min(lambda(kk),0.)/SL)-(max(lambda(kk),0.)/SR)
end do

do kk=1,7 
  if(SL>=0.)then
    flux(jj,kk)=fleft(kk)
  else if(SR<=0.)then
    flux(jj,kk)=fright(kk)
  else if(SL<0. .and. SR>0.)then
    flux(jj,kk)=SR*fleft(kk)-SL*fright(kk)+SL*SR*(uright(kk)-uleft(kk))
    !HLLI correction terms
    if(HLLIentropymode==1)then
    flux(jj,kk)=flux(jj,kk)-(SR*SL)*delta(4)*R4(kk)*alpha(4)
    end if
    if(HLLIslowmode==1)then
    flux(jj,kk)=flux(jj,kk)-(SR*SL)&
    *(delta(3)*R3(kk)*alpha(3)+delta(5)*R5(kk)*alpha(5)) 
    end if
    if(HLLIalfvenmode==1)then
    flux(jj,kk)=flux(jj,kk)-(SR*SL)&
    *(delta(2)*R2(kk)*alpha(2)+delta(6)*R6(kk)*alpha(6)) 
    end if
    if(HLLIfastmode==1)then
    flux(jj,kk)=flux(jj,kk)-(SR*SL)&
    *(delta(1)*R1(kk)*alpha(1)+delta(7)*R7(kk)*alpha(7)) 
    end if
    flux(jj,kk)=flux(jj,kk)/(SR-SL)
  end if
end do

!-------------------------------------------------------------
!Cell interface electric field
!-------------------------------------------------------------
if(dir==1)then
  fx(jj,ind)=flux(jj,5)
else if(dir==2)then
  fy(ind,jj)=flux(jj,5)
end if


if(debug==1)then
  print*,'flux=',flux(jj,1),flux(jj,2),flux(jj,3), &
                 flux(jj,4),flux(jj,5),flux(jj,6),flux(jj,7)
end if 

end do
!$omp end parallel do

end subroutine computeFluxHLLI


subroutine timeStep()
integer::ii,jj
real*8::att

!Compute max wave speed. Determine Roe averaged states 
!at each cell interface in both x and y-directions.
!Compute maximum wave propagation speed from among these.

!reset smax
!smax=0._8
!print*,'Start computing time-step.'

!x-sweeps
!$omp parallel do & 
!$omp& private(ii,rhoL,uL,vL,wL,bxL,byL,bzL,BBL,vvL,pL, &
!$omp& ptotL,fleft,rhoR,uR,vR,wR,bxR,byR,bzR,BBR,vvR,pR, & 
!$omp& ptotR,fright,uleft,uright,rhot,ut,vt,wt,pt,byt,bzt,at, &
!$omp& att,ptot,bx,by,bz,bsqr,cf) reduction(max:smax)
do jj=1,ny
  do ii=1,nx

   rhoL=u1(ii,jj,1)
   uL=u1(ii,jj,2)/u1(ii,jj,1)
   vL=u1(ii,jj,3)/u1(ii,jj,1)
   wL=u1(ii,jj,4)/u1(ii,jj,1)
   bxL=u1(ii,jj,5)
   byL=u1(ii,jj,6)
   bzL=u1(ii,jj,7)
   BBL=bxL*bxL+byL*byL+bzL*bzL
   vvL=uL*uL+vL*vL+wL*wL
   pL=(gam-1.)*(u1(ii,jj,8)-0.5_8*rhoL*vvL-0.5_8*BBL)
   ptotL=pL+0.5_8*BBL


   rhoR=u1(ii+1,jj,1)
   uR=u1(ii+1,jj,2)/u1(ii+1,jj,1)
   vR=u1(ii+1,jj,3)/u1(ii+1,jj,1)
   wR=u1(ii+1,jj,4)/u1(ii+1,jj,1)
   bxR=u1(ii+1,jj,5)
   byR=u1(ii+1,jj,6)
   bzR=u1(ii+1,jj,7)
   BBR=bxR*bxR+byR*byR+bzR*bzR
   vvR=uR*uR+vR*vR+wR*wR
   pR=(gam-1.)*(u1(ii+1,jj,8)-0.5_8*rhoR*vvR-0.5_8*BBR)
   ptotR=pR+0.5_8*BBR
 
   !print*,'i,j=',ii,jj
   !print*,'rhot,pt',rhot,pt
 

   !----------------------------------------------
   !Roe-Averaged State Variables (arithmetic mean):
   !----------------------------------------------
   rhot=0.5_8*(rhoL+rhoR)
   ut=0.5_8*(uL+uR)
   byt=0.5_8*(byL+byR)
   bzt=0.5_8*(bzL+bzR)
   ptot=0.5_8*(ptotL+ptotR) !total pressure: fluid+magnetic field
   pt=ptot-0.5_8*(bxL*bxL+byt*byt+bzt*bzt)
   
   if(pt>0._8 .and. rhot>0._8)then
   at=sqrt(gam*pt/rhot)


   !-------------------------------------------------
   !Fast Mode wave speedc_f
   !-------------------------------------------------
   bx=sqrt(bxL*bxL/rhot)
   by=sqrt(byt*byt/rhot)
   bz=sqrt(bzt*bzt/rhot)
   bsqr=(bxL*bxL+byt*byt+bzt*bzt)/rhot

   cf=0.5_8*(at*at+bsqr+sqrt((at*at+bsqr)*(at*at+bsqr)-4._8*at*at*bx*bx)) 
   cf=sqrt(cf)

   !find maximum speed
   smax=max(smax,max(abs(ut),abs(vt),abs(wt))+cf)

   end if

  
  end do 
end do
 !$omp end parallel do

!y-sweeps
!$omp parallel do private(jj,rhoL,uL,vL,wL,bxL,byL,bzL,BBL,vvL,pL, &
!$omp& ptotL,fleft,rhoR,uR,vR,wR,bxR,byR,bzR,BBR,vvR,pR, & 
!$omp& ptotR,fright,uleft,uright,rhot,ut,vt,wt,pt,byt,bzt,at, &
!$omp& att,ptot,bx,by,bz,bsqr,cf) reduction(max:smax)
do ii=1,nx
  do jj=1,ny
  
   rhoL=u1(ii,jj,1)
   uL=u1(ii,jj,3)/u1(ii,jj,1)
   vL=u1(ii,jj,2)/u1(ii,jj,1)
   wL=u1(ii,jj,4)/u1(ii,jj,1)
   bxL=u1(ii,jj,6)
   byL=u1(ii,jj,5)
   bzL=u1(ii,jj,7)
   BBL=bxL*bxL+byL*byL+bzL*bzL
   vvL=uL*uL+vL*vL+wL*wL
   pL=(gam-1.)*(u1(ii,jj,8)-0.5_8*rhoL*vvL-0.5*BBL)
   ptotL=pL+0.5_8*BBL

   rhoR=u1(ii,jj+1,1)
   uR=u1(ii,jj+1,3)/u1(ii,jj+1,1)
   vR=u1(ii,jj+1,2)/u1(ii,jj+1,1)
   wR=u1(ii,jj+1,4)/u1(ii,jj+1,1)
   bxR=u1(ii,jj+1,6)
   byR=u1(ii,jj+1,5)
   bzR=u1(ii,jj+1,7)
   BBR=bxR*bxR+byR*byR+bzR*bzR
   vvR=uR*uR+vR*vR+wR*wR
   pR=(gam-1.)*(u1(ii,jj+1,8)-0.5_8*rhoR*vvR-0.5*BBR)
   ptotR=pR+0.5_8*BBR
 
   !----------------------------------------------
   !Roe-Averaged State Variables (arithmetic mean):
   !----------------------------------------------
   rhot=0.5_8*(rhoL+rhoR)
   ut=0.5_8*(uL+uR)
   byt=0.5_8*(byL+byR)
   bzt=0.5_8*(bzL+bzR)
   ptot=0.5_8*(ptotL+ptotR) !total pressure: fluid+magnetic field
   pt=ptot-0.5_8*(bxL*bxL+byt*byt+bzt*bzt)
   
   if(pt>0._8 .and. rhot>0._8)then

   at=sqrt(gam*pt/rhot)

   
   !-------------------------------------------------
   !Fast Mode wave speedc_f
   !-------------------------------------------------
   bx=sqrt(bxL*bxL/rhot)
   by=sqrt(byt*byt/rhot)
   bz=sqrt(bzt*bzt/rhot)
   bsqr=(bxL*bxL+byt*byt+bzt*bzt)/rhot

   cf=0.5_8*(at*at+bsqr+sqrt((at*at+bsqr)*(at*at+bsqr)-4._8*at*at*bx*bx))
   cf=sqrt(cf)

   !find maximum speed
   smax=max(smax,max(abs(ut),abs(vt),abs(wt))+cf)
   end if

  end do
end do 
!$omp end parallel do
 
!--------------------------------------------
!Set time-step size
!--------------------------------------------
dt=dx*cour/smax 

!print*,'Done computing time-step.'

end subroutine timeStep


subroutine horcut()
!Local Variables
integer,parameter::ycell=ny/2
integer::jj

do jj=1,nx
  x=xmin+(jj-0.5)*dx
  dens=u2(jj,ycell,1)
  velx=u1(jj,ycell,2)/u1(jj,ycell,1)
  vely=u1(jj,ycell,3)/u1(jj,ycell,1)
  velz=u1(jj,ycell,4)/u1(jj,ycell,1)
  magx=u1(jj,ycell,5)
  magy=u1(jj,ycell,6)
  magz=u1(jj,ycell,7)
  pres=(gam-1.)*(u1(jj,ycell,8)-0.5*dens*(velx**2.+vely**2.+velz**2.)&
        -0.5*(magx**2.+magy**2.+magz**2.))   
  
  write(12,*) x,dens,velx,vely,pres
  !write(13,*) x,magx,magy,magz
end do

end subroutine horcut

subroutine vercut()
!Local Variables
integer,parameter::xcell=nx/2
integer::jj

do jj=1,nx
  y=ymin+(jj-0.5)*dy
  dens=u1(xcell,jj,1)
  velx=u1(xcell,jj,2)/u1(xcell,jj,1)
  vely=u1(xcell,jj,3)/u1(xcell,jj,1)
  velz=u1(xcell,jj,4)/u1(xcell,jj,1)
  magx=u1(xcell,jj,5)
  magy=u1(xcell,jj,6)
  magz=u1(xcell,jj,7)
  pres=(gam-1.)*(u1(xcell,jj,8)-0.5*dens*(velx**2.+vely**2.+velz**2.)&
        -0.5*(magx**2.+magy**2.+magz**2.))   
  
  write(14,*) y,dens,velx,vely,velz,pres
  write(15,*) y,magx,magy,magz
end do

end subroutine vercut

subroutine diagcut()
integer::jj
real*8::vperp,vpar,bpar,bperp

do jj=1,nx
  y=ymin+(jj-0.5)*dy
  dens=u1(jj,jj,1)
  velx=u1(jj,jj,2)/u1(jj,jj,1)
  vely=u1(jj,jj,3)/u1(jj,jj,1)
  velz=u1(jj,jj,4)/u1(jj,jj,1)
  vpar=(vely+velx)/sqrt(2.)
  vperp=(vely-velx)/sqrt(2.)
  magx=u1(jj,jj,5)
  magy=u1(jj,jj,6)
  magz=u1(jj,jj,7)
  bpar=(magy+magx)/sqrt(2.)
  bperp=(magy-magx)/sqrt(2.)
  pres=(gam-1.)*(u1(jj,jj,8)-0.5*dens*(velx**2.+vely**2.+velz**2.)&
        -0.5*(magx**2.+magy**2.+magz**2.))   
  
  write(16,*) y,dens,vpar,vperp,velz,pres
  write(17,*) y,bpar,bperp,magz
end do

end subroutine diagcut

subroutine fileOutput(iunit,cell_count)
!Input variables
integer::iunit
integer::cell_count(-1:nx+2,-1:ny+2)
!local variables
integer::jj,kk
real*8::magpres

do kk=1,ny
  y=ymin+(kk-0.5)*dy
  do jj=1,nx
   x=xmin+(jj-0.5)*dx
   dens=u2(jj,kk,1)
   velx=u2(jj,kk,2)/u2(jj,kk,1)
   vely=u2(jj,kk,3)/u2(jj,kk,1)
   velz=u2(jj,kk,4)/u2(jj,kk,1)
   magx=u2(jj,kk,5)
   magy=u2(jj,kk,6)
   magz=u2(jj,kk,7)
   magpres=0.5*(magx**2.+magy**2.+magz**2.)
   pres=(gam-1.)*(u2(jj,kk,8)-0.5*dens*(velx**2.+vely**2.+velz**2.)&
        -0.5*(magx**2.+magy**2.+magz**2.))    

   write(iunit,*) x,y,(dens),pres,vely!,real(cell_count(jj,kk)),pres,magpres
   !if(mod(jj,2)==0 .and. mod(kk,2)==0)then
   !   write(iunit+7000,*) x,y,magx/sqrt(magx**2.+magy**2.+magz**2.),&
   !   magy/sqrt(magx**2.+magy**2.+magz**2.)
   !end if
  end do
end do

end subroutine fileOutput


subroutine bound()
!Local Variables
integer::jj,kk,ll

if (boundaryType==1)then
!$omp parallel do shared(u2) private(ll)
do jj=-1,nx+2
  do ll=1,8
    u2(jj,0,ll)=u2(jj,1,ll)	
    u2(jj,-1,ll)=u2(jj,1,ll)
    u2(jj,ny+1,ll)=u2(jj,ny,ll)
    u2(jj,ny+2,ll)=u2(jj,ny,ll)      
  end do
end do
!$omp end parallel do

!$omp parallel do shared(u2) private(ll)
do kk=-1,ny+2
  do ll=1,8
    u2(0,kk,ll)=u2(1,kk,ll)	
    u2(-1,kk,ll)=u2(1,kk,ll)
    u2(nx+1,kk,ll)=u2(nx,kk,ll)
    u2(nx+2,kk,ll)=u2(nx,kk,ll)
  end do
end do
!$omp end parallel do

else if(boundaryType==2)then
!$omp parallel do shared(u2) private(ll)
do jj=-1,nx+2
  do ll=1,8
    u2(jj,0,ll)=u2(jj,1,ll)	
    u2(jj,-1,ll)=u2(jj,1,ll)
    u2(jj,ny+1,ll)=u2(jj,ny,ll)
    u2(jj,ny+2,ll)=u2(jj,ny,ll)      
  end do
end do
!$omp end parallel do

!$omp parallel do shared(u2) private(ll)
do kk=-1,ny+2
  do ll=1,8
    u2(0,kk,ll)=u2(nx,kk,ll)	
    u2(-1,kk,ll)=u2(nx-1,kk,ll)
    u2(nx+1,kk,ll)=u2(1,kk,ll)
    u2(nx+2,kk,ll)=u2(2,kk,ll)
  end do
end do
!$omp end parallel do

else if(boundaryType==3)then
!$omp parallel do shared(u2) private(ll)
do jj=-1,nx+2
  do ll=1,8
    u2(jj,0,ll)=u2(jj,1,ll)	
    u2(jj,-1,ll)=u2(jj,1,ll)
    u2(jj,ny+1,ll)=u2(jj,ny,ll)
    u2(jj,ny+2,ll)=u2(jj,ny,ll)      
  end do
  u2(jj,0,3)=-u2(jj,1,3)	
  u2(jj,-1,3)=-u2(jj,1,3)
  u2(jj,ny+1,3)=-u2(jj,ny,3)
  u2(jj,ny+2,3)=-u2(jj,ny,3)
  u2(jj,0,6)=-u2(jj,1,6)	
  u2(jj,-1,6)=-u2(jj,1,6)
  u2(jj,ny+1,6)=-u2(jj,ny,6)
  u2(jj,ny+2,6)=-u2(jj,ny,6)

end do
!$omp end parallel do

!$omp parallel do shared(u2) private(ll)
do kk=-1,ny+2
  do ll=1,8
    u2(0,kk,ll)=u2(1,kk,ll)	
    u2(-1,kk,ll)=u2(1,kk,ll)
    u2(nx+1,kk,ll)=u2(nx,kk,ll)
    u2(nx+2,kk,ll)=u2(nx,kk,ll)
  end do
  u2(0,kk,2)=-u2(1,kk,2)	
  u2(-1,kk,2)=-u2(1,kk,2)
  u2(nx+1,kk,2)=-u2(nx,kk,2)
  u2(nx+2,kk,2)=-u2(nx,kk,2)
  u2(0,kk,5)=-u2(1,kk,5)	
  u2(-1,kk,5)=-u2(1,kk,5)
  u2(nx+1,kk,5)=-u2(nx,kk,5)
  u2(nx+2,kk,5)=-u2(nx,kk,5)
end do
!$omp end parallel do

else if(boundaryType==4)then
!$omp parallel do shared(u2) private(ll)
do jj=-1,nx+2
  do ll=1,8
    u2(jj,0,ll)=u2(jj,ny,ll)	
    u2(jj,-1,ll)=u2(jj,ny-1,ll)
    u2(jj,ny+1,ll)=u2(jj,1,ll)
    u2(jj,ny+2,ll)=u2(jj,2,ll)      
  end do
end do
!$omp end parallel do

!$omp parallel do shared(u2) private(ll)
do kk=-1,ny+2
  do ll=1,8
    u2(0,kk,ll)=u2(nx,kk,ll)	
    u2(-1,kk,ll)=u2(nx-1,kk,ll)
    u2(nx+1,kk,ll)=u2(1,kk,ll)
    u2(nx+2,kk,ll)=u2(2,kk,ll)
  end do
end do
!$omp end parallel do

else if(boundaryType==5)then
!$omp parallel do shared(u2) private(ll)
do jj=-1,nx+2
  do ll=1,8
    u2(jj,0,ll)=u2(jj,1,ll)	
    u2(jj,-1,ll)=u2(jj,1,ll)
    u2(jj,ny+1,ll)=u2(jj,ny,ll)
    u2(jj,ny+2,ll)=u2(jj,ny,ll)      
  end do
  u2(jj,0,3)=-u2(jj,1,3)	
  u2(jj,-1,3)=-u2(jj,1,3)
  u2(jj,ny+1,3)=-u2(jj,ny,3)
  u2(jj,ny+2,3)=-u2(jj,ny,3)
  u2(jj,0,6)=-u2(jj,1,6)	
  u2(jj,-1,6)=-u2(jj,1,6)
  u2(jj,ny+1,6)=-u2(jj,ny,6)
  u2(jj,ny+2,6)=-u2(jj,ny,6)
end do
!$omp end parallel do

!$omp parallel do shared(u2) private(ll)
do kk=-1,ny+2
  do ll=1,8
    u2(0,kk,ll)=u2(nx,kk,ll)	
    u2(-1,kk,ll)=u2(nx-1,kk,ll)
    u2(nx+1,kk,ll)=u2(1,kk,ll)
    u2(nx+2,kk,ll)=u2(2,kk,ll)
  end do
end do
!$omp end parallel do

end if

end subroutine bound


subroutine protection()
!Local Variables
integer::jj,kk,ll
real*8::pres0

do jj=1,nx
  do kk=1,ny
   
    dens=u2(jj,kk,1)
    if(dens<0._8)then
      !if(debug==1)then
        print*,'Negative density detected in cell:',jj,kk,'. Using protection.'
      !end if
      dens=max(u1(jj,kk,1),min_dens)
      u2(jj,kk,1)=dens
    end if

    pres=(gam-1.)*( u2(jj,kk,8)- &
         0.5*(u2(jj,kk,2)**2.+u2(jj,kk,3)**2.+u2(jj,kk,4)**2.)/u2(jj,kk,1) &
        -0.5*(u2(jj,kk,5)**2.+u2(jj,kk,6)**2.+u2(jj,kk,7)**2.) )    
    pres0=(gam-1.)*( u1(jj,kk,8)- &
         0.5*(u1(jj,kk,2)**2.+u1(jj,kk,3)**2.+u1(jj,kk,4)**2.)/u1(jj,kk,1) &
        -0.5*(u1(jj,kk,5)**2.+u1(jj,kk,6)**2.+u1(jj,kk,7)**2.) ) 

    if(pres<0._8)then
    !if(debug==1)then
      print*,'Negative pressure detected in cell:',jj,kk,'. Using protection.'
      print*,'pres=',pres
     ! end if
      pres=max(min_pres,pres0)
      print*,'new value=',pres
      u2(jj,kk,8)=gamul*pres &
        +0.5*(u1(jj,kk,2)**2.+u1(jj,kk,3)**2.+u1(jj,kk,4)**2.)/u1(jj,kk,1) &
        +0.5*(u1(jj,kk,5)**2.+u1(jj,kk,6)**2.+u1(jj,kk,7)**2.) 
    end if
   end do
end do

end subroutine protection

function sgn(x) result(fx)
  real*8,intent(in)::x
  real*8::fx
  fx=sign(1._8,x)
end

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


end module RoeSolver2D_mod
