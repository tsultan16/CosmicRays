!Chang-Cooper(1970) algorithm for Isotropic Cosmic Ray-Fokker Planck Equation

module FokkerPlanck_mod
use constants_mod
implicit none


contains

subroutine CRinit()

integer::i,j

open(unit=10,file='fp.txt')
open(unit=13,file='fx.txt')
!----------------------------------------------
!Initialization
!----------------------------------------------
!Momentum space bounds
pmin=1.D-3
pmax=1.D3

print*,'np,pmin,pmax=',np,pmin,pmax

nparticle=0.

!Generate momentum mesh points
do i=-2,np+3
  if(meshType==1)then
   px(i)=pmin+i*(pmax-pmin)/np
  else if(meshType==2)then
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
  f(i,j)=1.d-12*px(i)**(-5.) !exp(-( (xmin+j*dx-0.5)/(10.*dx) )**2. )
 end do
 f(-1,j)=f(0,j)  
 f(np,j)=f(np-1,j)
end do


fp=0.
do i=0,np-1
 do j=0,nx-1
  fp(i)=fp(i)+dx*f(i,j)
 end do
 write(10,*) px(i),f(i,5)!+1.D-80
 !print*,'p,f(p)',px(i),fp(i)
end do

fx=0.
do j=0,nx-1
  do i=0,np-1
  fx(j)=fx(j)+dp(i)*f(i,j)
 end do
 write(13,*) xmin+j*dx,f(50,j)!fx(j)!+1.D-50
 !print*,'x,f(x)',xmin+j*dx,fx(j)
end do


print*,'Done initializing CR distribution.'

end subroutine

subroutine CRevolve()
!Evolve Distribution Function
integer::j,k

!-------------------------------------
!Mometum Operator Sweeps   
!-------------------------------------
  do j=0,nx-1
   ftemp1=0.
   do k=-1,np
    ftemp1(k)=f(k,j)
   end do
    
    if(methodType==3)then
      call computePFlux(ftemp1,j)  

    else
      call computeMomentumCoefficients(j)
      call tridiagSolveMomentum(np,ftemp1)
      !print*,'CHECKPOINT1'
      !call momentumAdvect(j,ftemp1)
      !print*,'CHECKPOINT2'
    end if

   ftemp1(-1)=ftemp1(0)  
   ftemp1(np)=ftemp1(np-1)

   do k=-1,np
    f(k,j)=ftemp1(k)  
    !print*,'j,k,ftemp1=',j,k,ftemp1(k) 
   end do

  end do

!-------------------------------------
!Spatial Operator Sweeps   
!------------------------------------- 
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


subroutine momentumAdvect(xk,ftemp)
!Input Variables
integer::xk
real*8::ftemp(-1:np)

!Local Variables
integer::jj,pflux(-1:np+2)
real*8::B2,temp

!-------------------------------------------------------------
! F^n_i+1/2 (first order upwind flux)
!-------------------------------------------------------------
do jj=-1,np-1
 B2=B(jj,xk)
 !print*,'jj,B_jj+1/2,ftemp(jj)=',jj,B2,ftemp(jj)
 if(B2>0._8)then
  pflux(jj)=B2*ftemp(jj+1)
 else
  temp=B2*ftemp(jj)
  print*,'CHECKPOINT3,jj,B2,ftemp(jj),temp=',jj,B2,ftemp(jj),temp
  pflux(jj)=temp
  print*,'CHECKPOINT4'
 end if
 !print*,'jj,pflux(jj)=',jj,pflux(jj)
end do

!No-flux Boundary Condition
pflux(-1)=0.
pflux(np-1)=0.

!Update distribution function
!fupdate=ftemp1
do jj=0,np-1
ftemp(jj)=ftemp(jj)*(1.+dt*Hp(jj)) &
           +(dt/A(jj))*(pflux(jj)-pflux(jj-1))/dp(jj) &
           + dt*S(jj,t)
end do

end subroutine momentumAdvect



subroutine computeMomentumCoefficients(xk)

!Input variables
integer::xk

integer::jj
real*8::w(-2:np+2),Wplus(-1:np+1),Wminus(-1:np+1),Wbig(-1:np+1)
real*8::B2,C2,B1,C1,pi2,pi2L,w1,w2

!-------------------------------------------------------------
!Weights: w_j, W_j, W^+_j+1/2, W^-_j+1/2 
!-------------------------------------------------------------
if(methodType .ne. 4)then
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
  !print*,'j,w,Wbig=',j,w(j),Wbig(j)
end do

do jj=-1,np+1
  Wminus(jj)=0.5*(Wbig(jj)*exp(-0.5*w(jj))+Wbig(jj+1)*exp(-0.5*w(jj+1)) )
  Wplus(jj)=0.5*(Wbig(jj)*exp(0.5*w(jj))+Wbig(jj+1)*exp(0.5*w(jj+1)) )
end do


!No-flux Boundary condition requirement
Wplus(-1)=0.
Wminus(-1)=0.
Wplus(np-1)=0.
Wminus(np-1)=0.

!At=0.
!Bt=0.
!Ct=0.

end if

!------------------------------------------------------------
!Co-efficients for tridiagonal system: ~A_j,~B_j,~C_j
!------------------------------------------------------------
do jj=0,np-1
  C1=0.5*(C(jj-1)+C(jj))  !C_j-1/2
  C2=0.5*(C(jj)+C(jj+1))  !C_j+1/2
  B1=0.5*(B(jj-1,xk)+B(jj,xk))
  B2=0.5*(B(jj,xk)+B(jj+1,xk))
  w1=0.5*(w(jj-1)+w(jj))
  w2=0.5*(w(jj)+w(jj+1))

  if(methodType .ne. 4)then
   At(jj)=dt/(dp(jj)*dp2(jj-1)*A(jj))
   At(jj)=At(jj)*C1*Wminus(jj-1)
   At(jj)=-At(jj)
  end if
  if(methodType==2)then
   At(jj)=At(jj)*0.5
  end if

  if(methodType==4)then
  At(jj)=dt/(dp(jj)*dp2(jj-1)*A(jj))
  At(jj)=At(jj)*(C1-0.5*B1*dp2(jj-1))
  At(jj)=-At(jj)
  end if

  if(methodType .ne. 4)then
   Bt(jj)= C1*Wplus(jj-1)/dp2(jj-1) &
        +C2*Wminus(jj)/dp2(jj)
   Bt(jj)=Bt(jj)*dt/(dp(jj)*A(jj))
   Bt(jj)=Bt(jj)+1.+dt*Hp(jj)
  end if
  if(methodType==2)then 
   Bt(jj)=0.5*(Bt(jj)+1.)
  end if

  if(methodType==4)then
   Bt(jj)=(C1+0.5*B1*dp2(jj-1))/dp2(jj-1) &
        +(C2+0.5*B2*dp2(jj))/dp2(jj)
   Bt(jj)=Bt(jj)*dt/(dp(jj)*A(jj))
   Bt(jj)=Bt(jj)+1.+dt*Hp(jj)
  end if 

  if(methodType .ne. 4)then
   Ct(jj)=dt/(dp(jj)*dp2(jj)*A(jj))
   Ct(jj)=Ct(jj)*C2*Wplus(jj)
   Ct(jj)=-Ct(jj)
  end if
  if(methodType==2)then
   Ct(jj)=Ct(jj)*(1.+0.5*w2)
  end if
  if(methodType==4)then
   Ct(jj)=dt/(dp(jj)*dp2(jj)*A(jj))
   Ct(jj)=Ct(jj)*(C2+0.5*B2*dp2(jj))
   Ct(jj)=-Ct(jj)
  end if

  !print*,'jj,At,Bt,Ct=',jj,At(jj),Bt(jj),Ct(jj)
end do
 
end subroutine computeMomentumCoefficients

subroutine tridiagSolveMomentum(n,ftemp) !from Numerical Recipes book

!Input variables
integer::n
real*8::ftemp(-1:n)

real*8::bet,gam(0:n-1),r(0:n-1)

integer::kk


do kk=0,n-1
 !print*,'k,p,S(p,t)=',kk,px(kk),S(kk,t)
 r(kk)=ftemp(kk)+dt*S(kk,t)
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
real*8::w(-2:nx+2),Wplus(-1:nx+1),Wminus(-1:nx+1),Wbig(-1:nx+1)
real*8::B2,C2,B1,C1,pi2,pi2L

!-------------------------------------------------------------
!Weights: w_j, W_j, W^+_j+1/2, W^-_j+1/2 
!-------------------------------------------------------------
do jj=-1,nx+2
 !print*,'jj=',jj
 B2=E(jj)
 C2=Kxx(pk)
 w(jj)=dx*B2/C2 

 go to 123
 if(abs(w(jj))>1.d4 .or. abs(w(jj))<1.d-12)then
   print*,'***ERROR***, w_j=',w(jj)
   print*,'dp,B2,C2=',dp(pk),B2,C2
   print*,'w_j is out of bounds.' 
   print*,'Need to re adjust size of momentum advection and diffusion co-efficients.' 
   STOP
 end if
 123 continue

  Wbig(jj)=0.5*w(jj)/sinh(0.5*w(jj))
  !print*,'jj,w,Wbig=',jj,w(jj),Wbig(jj)
end do

do jj=-1,nx+1
   Wminus(jj)=0.5*(Wbig(jj)*exp(-0.5*w(jj))+Wbig(jj+1)*exp(-0.5*w(jj+1)) )
   Wplus(jj)=0.5*(Wbig(jj)*exp(0.5*w(jj))+Wbig(jj+1)*exp(0.5*w(jj+1)) )
end do


!No-flux Boundary condition requirement
Wplus(-1)=0.
Wminus(-1)=0.
Wplus(nx-1)=0.
Wminus(nx-1)=0.


!------------------------------------------------------------
!Co-efficients for tridiagonal system: ~A_j,~B_j,~C_j
!------------------------------------------------------------
do jj=0,nx-1
  C1=Kxx(pk)
  C2=C1

  Atx(jj)=dt/(dx*dx)
  Atx(jj)=Atx(jj)*C1*Wminus(jj-1)
  Atx(jj)=-Atx(jj)
  if(methodType==2)then
   Atx(jj)=Atx(jj)*0.5 
  end if
  

  Btx(jj)= C1*Wplus(jj-1)/dx &
        +C2*Wminus(jj)/dx
  Btx(jj)=Btx(jj)*dt/dx
  Btx(jj)=Btx(jj)+1.+dt*Hx(pk)
  if(methodType==2)then 
   Btx(jj)=0.5*(Btx(jj)+1.)
  end if

  Ctx(jj)=dt/(dx*dx)
  Ctx(jj)=Ctx(jj)*C2*Wplus(jj)
  Ctx(jj)=-Ctx(jj)
  if(methodType==2)then
   Ctx(jj)=Ctx(jj)*0.5
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
 r(kk)=ftemp(kk)+dt*L(kk,t)
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

 !bx=0.00001*(px(pk)**3.)*dudx(xk)/3.
  bx=(px(pk)**3.)*dudx(xk)/3.


 !if(abs(bx)<1.d-10) bx=1.d-10 

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
function Kxx(pk) result(kx)

integer,intent(in)::pk
real*8::kx

kx=1.*(px(pk)**0.6)

end function Kxx


!Momentum Source term co-efficient
function Hp(xk) result(hx)
 
integer,intent(in)::xk
real*8::hx

hx=0.!-dudx/3.

end function Hp


!Spatial Source term co-efficient
function Hx(pk) result(hxx)
 
integer,intent(in)::pk
real*8::hxx

hxx=0.!dudx

end function Hx


!Momentum Space Injection term
function S(pk,t) result(sx)
 
integer,intent(in)::pk
real*8,intent(in)::t
real*8::sx

sx=0.

end function S

!Spatial Injection term
function L(xk,t) result(lx)
 
integer,intent(in)::xk
real*8,intent(in)::t
real*8::lx

!lx=exp(-( (xmin+xk*dx-0.8)/(dx) )**2. )
if(xk==shockCell)then
lx=0.
end if

end function L

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
           + dt*S(jj,t)
end do
!ftemp1=fupdate

end subroutine computePFlux


end module FokkerPlanck_mod
