!Chang-Cooper(1970) algorithm for Cosmic Ray-Fokker Planck Equation
!(evolution in momentum space only+isotropic plasma)


program FokkerPlanck_Implicit

implicit none

integer,parameter::np=100 !momentum mesh size
integer,parameter::nx=10 !spatial mesh size
integer,parameter::nt=200 !time steps
integer,parameter::meshType=2 !1:uniform, 2:logarithmic
integer,parameter::testType=0 !0,1,2,3


real,parameter::tol=1.D-8
real*8::f(-1:np+2,-1:nx+2),ftemp1(0:np-1),ftemp2(0:nx)
real*8::At(0:np-1),Bt(0:np-1),Ct(0:np-1)
real*8::px(-2:np+3),dp(-1:np+1),dp2(-1:np+1)
real*8::pmin,pmax,t,nparticle
real*8::xmin,xmax,dx
real*8::dt=0.001
!-----------------------------
real*8::ek(0:np-1),gk(0:np-1)
!-----------------------------
integer::i,j,k

open(unit=10,file='fp.txt')
print*,'Simple Fully Implicit Fokker-Planck Solver.'

!----------------------------------------------
!Initialization
!----------------------------------------------
!Momentum space bounds
pmin=1.D-3
pmax=1.D3
xmin=0.
xmax=1.

dx=(xmax-xmin)/nx


print*,'np,pmin,pmax=',np,pmin,pmax


!Generate momentum mesh points
do i=-2,np+2
  if(meshType==1)then
   px(i)=pmin+i*(pmax-pmin)/np
  else if(meshType==2)then
   px(i)=pmin*(pmax/pmin)**(real(i)/real(np))
  end if
end do

!Momentum mesh spacing
do i=-1,np+1
 dp(i)=0.5*(px(i+1)-px(i-1))
 dp2(i)=px(i+1)-px(i)
end do


 
do i=0,np-1
  do j=-1,nx+1
  !Initial profile
  !if(testType==3)then 
  !f(i,j)=exp(-((px(i)-(pmax-pmin)/2.)/(5.*dp(i)))**2.)
  
  !else if(testType==0)then
  !  f(i,j)=1.D-15*px(i)**(-5.) 
  !else
    f(i,j)=0.!+1.D-60 
  !end if

end do

!f(-1,j)=0.
!f(np+1,j)=0.
end do

nparticle=0.
do i=0,np-1
 nparticle=nparticle+f(i,1)
end do

print*,'nparticle=',nparticle

do i=0,np-1
  write(10,*) px(i),f(i,1)+1.D-70
end do

!----------------------------------------
!Evolve Distribution Function
!----------------------------------------
t=0._8
do i=1,nt
  print*,'Time Step=',i

  !-------------------------------------
  !Mometum Operator Sweeps   
  !-------------------------------------
  do j=0,nx

   do k=0,np-1
    ftemp1(k)=f(k,j)
   end do

   call computeMomentumCoefficients(j)
   call tridiagSolve(np,ftemp1)
   !call thomasAlgorithm(np,ftemp1)   

   do k=0,np-1
    f(k,j)=ftemp1(k)
   end do

  end do

  !-------------------------------------
  !Spatial Operator Sweeps   
  !------------------------------------- 
  do j=0,np

  end do


  do j=0,np-1
   !print*,'p,f(p)=',px(j),u(j)
    write(10,*) px(j),f(j,1) 
  end do


  t=t+i*dt
  print*,'% Complete=',(real(i)/real(nt))*100.
  !if(i<0.1*nt)dt=1.2*dt

 nparticle=0.
 do j=0,np-1
  nparticle=nparticle+f(j,1)
 end do

print*,'nparticle=',nparticle

end do

close(unit=10)
print*,'DONE.'

contains


subroutine computeMomentumCoefficients(xindex)

!Input variables
integer::xindex

integer::j
real*8::w(-1:np+2),Wplus(-1:np+2),Wminus(-1:np+2),Wbig(-1:np+2)
real*8::B2,C2,B1,C1,pi2,pi2L

!-------------------------------------------------------------
!Weights: w_j+1/2, W_j, W^+_j+1/2, W^-_j+1/2 
!-------------------------------------------------------------
do j=-1,np
 B2=0.5*( B(px(j),0._8,0._8)+B(px(j+1),0._8,0._8) )
 C2=0.5*( C(px(j),0._8)+C(px(j+1),0._8) )
 !print*,'j,dp2,B2,C2=',j,dp2(j),B2,C2
 w(j)=dp2(j)*B2/C2
 
end do

!------------------------------------------------------------
!Co-efficients for tridiagonal system: ~A_j,~B_j,~C_j
!------------------------------------------------------------
do j=0,np-1
  C1=0.5*( C(px(j-1),0._8)+C(px(j),0._8) )  !C_j-1/2
  C2=0.5*( C(px(j),0._8)+C(px(j+1),0._8) )  !C_j+1/2

  At(j)=dt/(dp(j)*dp2(j-1)*A(px(j)))
  At(j)=At(j)*C1*(1.-0.5*w(j-1))
  At(j)=-At(j)
  
  Bt(j)= C1*((1.+0.5*w(j-1)))/dp2(j-1) &
        +C2*(1.-0.5*w(j+1))/dp2(j)
  Bt(j)=Bt(j)*dt/(dp(j)*A(px(j)))
  Bt(j)=Bt(j)+1.-dt*Hp(px(j),0._8,0._8)

  Ct(j)=dt/(dp(j)*dp2(j)*A(px(j)))
  Ct(j)=Ct(j)*C2*(1.+0.5*w(j+1))
  Ct(j)=-Ct(j)
  
end do

!Impose No-flux boundary conditions
At(0)=0.
Ct(np-1)=0. 


!------------------------------------------------------------
!Compute Thomas algorithm co-efficients: e_k,g_k
!------------------------------------------------------------ 
!e(0)=0.
!g(0)=0. !These correspond to the boundary condition f(0)=0.

do j=0,np-1
 ek(j)=-Ct(j)/(Bt(j)+ek(j-1)*At(j))
 gk(j)=(f(j,xindex)-At(j)*gk(j-1))/(Bt(j)+ek(j-1)*At(j))
end do



end subroutine computeMomentumCoefficients

subroutine thomasAlgorithm(n,ftemp)

!Input variables
integer::n
real*8::ftemp(0:n-1)

integer::k

!Back substitution
do k=n-1,0,-1
 ftemp(k)=ftemp(k+1)*ek(k)+gk(k)
end do




end subroutine thomasAlgorithm


subroutine tridiagSolve(n,ftemp) !from Numerical Recipes book

!Input variables
integer::n
real*8::ftemp(0:n-1)

real*8::bet,gam(0:n-1),r(0:n-1)
integer::k

do k=0,np-1
 r(k)=ftemp(k)+dt*S(k,t)
end do

bet=Bt(0)
ftemp(0)=r(0)/bet

!Upper triangular decomposition and forward substitution
do k=1,n-1
  gam(k)=Ct(k-1)/bet
  bet=Bt(k)-At(k)*gam(k)
  ftemp(k)=(r(k)-At(k)*ftemp(k-1))/bet
end do

!Back substitution
do k=n-2,0,-1
 ftemp(k)=ftemp(k)-gam(k+1)*ftemp(k+1)
end do

end subroutine tridiagSolve

!Momentum Phase factor
function A(p) result(ax)
real*8,intent(in)::p
real*8::ax
  
 ax=1.!p*p

end function A

!Momentum Advection co-efficient
function B(p,t,dudx) result(bx)
real*8,intent(in)::p,t,dudx
real*8::bx

 bx=-(1.+p)!-(p**3.)*dudx/3.

end function B

!Momentum Diffusion co-efficient
function C(p,t) result(cx)

real*8,intent(in)::p,t
real*8::cx,Dpp

 Dpp=1.
 cx=Dpp*p*p

end function C

!Spatial Phase co-efficient
function D(p,t) result(dx)

real*8,intent(in)::p,t
real*8::dx

dx=1.

end function D


!Spatial Advection co-efficient
function E(u,t) result(ex)

real*8,intent(in)::u,t
real*8::ex

ex=-u

end function E

!Spatial Diffusion co-efficient
function Kxx(p,t) result(kx)

real*8,intent(in)::p,t
real*8::kx

kx=0.1*(p**0.6)

end function Kxx


!Momentum Source term co-efficient
function Hp(p,dudx,t) result(hx)
 
real*8,intent(in)::p,dudx,t
real*8::hx

hx=-1.!-dudx/3.

end function Hp


!Spatial Source term co-efficient
function Hx(dudx,t) result(hxx)
 
real*8,intent(in)::dudx,t
real*8::hxx

hxx=dudx

end function Hx


!Injection term
function S(pk,t) result(sx)
 
integer,intent(in)::pk
real*8,intent(in)::t
real*8::sx

sx=exp(-((px(pk)-0.1)/(0.05*dp(pk)))**2.)

end function S


end program FokkerPlanck_Implicit
