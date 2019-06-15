!Chang-Cooper(1970) algorithm for Cosmic Ray-Fokker Planck Equation
!(evolution in momentum space only+isotropic plasma)


program FokkerPlanck_CC

implicit none

integer,parameter::np=500 !momentum bins
integer,parameter::nt=100 !time steps
integer,parameter::meshType=2 !1:uniform, 2:logarithmic
integer,parameter::testType=2 !0,1,2,3
real,parameter::tol=1.D-8
real*8::f(0:np+1),r(0:np+1)
real*8::delta(-1:np),At(0:np),Bt(0:np),Ct(0:np)
real*8::px(-1:np+2),dp(-1:np+2),dp2(-1:np+2)
!real*8::e(0:np-1),g(0:np-1)
real*8::pmin,pmax,t
real*8::dt=0.1

integer::i,j

open(unit=10,file='fp.txt')


!----------------------------------------------
!Initialization
!----------------------------------------------
!Momentum space bounds
pmin=1.D-3
pmax=1.D3

print*,'np,pmin,pmax=',np,pmin,pmax

!Mesh points
do i=-1,np+2
  if(meshType==1)then
   px(i)=pmin+i*(pmax-pmin)/np
  else if(meshType==2)then
   px(i)=pmin*(pmax/pmin)**(real(i)/real(np))
  end if
end do

!Mesh spacing
do i=0,np+1
 dp(i)=0.5*(px(i+1)-px(i-1))
 dp2(i)=px(i+1)-px(i)
end do


do i=1,np
  !Gaussian Initial Profile
  !f(i)=exp(-( (px(i)-0.5*(pmax-pmin))/(0.02*np*dp(np/2)) )**2.)
  !Square Initial Profile
  !if(i>0.4*np .and. i< 0.6*np)then
  !  f(i)=1.
  !else
  !  f(i)=0.
  !end if

  !Powerlaw Initial profile
  if(testType==0)then
   if(i>0.1*np .and. i<0.9*np)then
    f(i)=1000.*(px(i))**(-5.)
   else
    f(i)=0.
   end if
  else
   f(i)=0.0
  end if

  write(10,*) px(i),f(i)

end do

!No-particle boundary conditions
!f(1)=0.!f(1)
!f(np)=0.!f(np)


!----------------------------------------
!Evolve Distribution Function
!----------------------------------------
t=0.
do i=1,nt
  print*,'Time Step=',i
  call computeCoefficients()
 
  !print*,'CHECKPOINT'
  !f(np)=0.
  !do j=np-1,2,-1
  ! f(j)=f(j+1)*e(j)+g(j)
  !end do 

  call tridiagSolve()

  f(0)=0.!f(1)
  f(np+1)=0.!f(np)
  
  do j=1,np
   write(10,*) px(j),f(j)
  end do

  t=t+i*dt
  print*,'% Complete=',(real(i)/real(nt))*100.

end do

close(unit=10)

contains


subroutine computeCoefficients()
integer::k
real*8::w(np),Wj(np),pi2,pi2L

!--------------------------------------
!First compute weights: w_j, W_j
!--------------------------------------
do k=0,np
 pi2=0.5*(px(k)+px(k+1))
 w(k)=dp(k)*B(pi2,0._8)/C(pi2,0._8)
 if(w(k)<tol)then
  Wj(k)=1.
 else
  Wj(k)=w(k)/(exp(w(k))-1.)
 end if
end do



!do k=-1,np
! wj=dp*B((k+0.5)*dp,0._8)/C((k+0.5)*dp,0._8)
 !print*,'j,w=',k,wj
! if(abs(wj)<tol)then
!   delta(k)=0.5
! else
!  delta(k)=(1./wj)-(1./(exp(wj)-1.))
! end if
!end do

!print*,'CHECKPOINT1'
!------------------------------------------------------------
!Compute co-efficients for tridiagonal system: ~A_j,~B_j,~C_j
!------------------------------------------------------------
do k=1,np
!wj=dp*B((k+0.5)*dp,0._8)/C((k+0.5)*dp,0._8)
!At(k)=dt*C((k+0.5)*dp,0._8)/(dp*dp*A(k*dp))
!At(k)=At(k)*((1.-delta(k))*wj+1.)
!At(k)=-At(k)
  pi2=0.5*(px(k)+px(k+1))
  pi2L=0.5*(px(k-1)+px(k))  

  !print*,'k,pi2L,pi2=',k,pi2,pi2L

  At(k)=dt/(dp(k)*dp2(k)*A(px(k)))
  At(k)=At(k)*C(pi2,0._8)*Wj(k)*exp(w(k))
  At(k)=-At(k)


!Bt(k)=dt/(dp*A(k*dp))
!Bt(k)=Bt(k)*( (C((k+0.5)*dp,0._8)+C((k-0.5)*dp,0._8))/dp &
  !    +(1.-delta(k-1))*B((k-0.5)*dp,0._8)-delta(k)*B((k+0.5)*dp,0._8) )+1.

  Bt(k)= C(pi2,0._8)*Wj(k) &
        +C(pi2L,0._8)*Wj(k-1)*exp(w(k-1)) 
  Bt(k)=Bt(k)*dt/(dp(k)*dp2(k)*A(px(k)))+1.
  Bt(k)=Bt(k)+dt*D(px(k),0._8)


!Ct(k)=dt/(dp*A(k*dp))
!Ct(k)=Ct(k)*( C((k-0.5)*dp,0._8)/dp-delta(k-1)*B((k-0.5)*dp,0._8)   )
!Ct(k)=-Ct(k)

  Ct(k)=dt/(dp(k)*dp2(k)*A(px(k)))
  Ct(k)=Ct(k)*C(pi2L,0._8)*Wj(k-1)
  Ct(k)=-Ct(k)
end do
!print*,'CHECKPOINT2'
 
!------------------------------------------------------------
!Compute Thomas algorithm co-efficients: e_k,g_k
!------------------------------------------------------------ 
!e(0)=0.
!g(0)=0. !These correspond to the boundary condition f(0)=0.

!do k=1,np-1
! e(k)=At(k)/(Bt(k)-e(k-1)*Ct(k))
! g(k)=(f(k)+Ct(k)*g(k-1))/(Bt(k)-e(k-1)*Ct(k))
!end do

end subroutine computeCoefficients

subroutine tridiagSolve() !from Numerical Recipes book

real*8::bet,gam(1:np)
integer::k

do k=1,np
 r(k)=f(k)+dt*E(k,t)
end do

bet=Bt(1)
f(1)=r(1)/bet

!Upper triangular decomposition and forward substitution
do k=2,np
  gam(k)=At(k-1)/bet
  bet=Bt(k)-Ct(k)*gam(k)
  f(k)=(r(k)-Ct(k)*f(k-1))/bet
end do

!Back substitution
do k=np-1,1,-1
  f(k)=f(k)-gam(k+1)*f(k+1)
end do


end subroutine tridiagSolve

!Phase factor
function A(p) result(ax)
real*8,intent(in)::p
real*8::ax
  
if(testType==0)then
 ax=0.01*p*p
else if(testType==1 .or. testType==2 .or. testType==3)then
 ax=1.
end if



end function A

!Advection co-efficient
function B(p,t) result(bx)
real*8,intent(in)::p,t
real*8::bx

if(testType==0)then
 bx=-0.1 
else if(testType==1)then
 bx=-(1.+p)
else if(testType==2)then
 bx=-p
else if(testType==3)then
 bx=-p*p
end if

end function B

!Diffusion co-efficient
function C(p,t) result(cx)

real*8,intent(in)::p,t
real*8::cx

if(testType==0)then
 cx=10. 
else if(testType==1)then
 cx=p*p
else if(testType==2)then
 cx=p*p
else if(testType==3)then
 cx=p*p*p
end if

end function C

!Escape term co-efficient
function D(p,t) result(dx)

real*8,intent(in)::p,t
real*8::dx

if(testType==0)then
 dx=0. 
else if(testType==1)then
 dx=-1.
else if(testType==2)then
 dx=-1./p
else if(testType==3)then
 dx=-1.
end if

end function D


!Source/Injection term co-efficient
function E(pk,t) result(ex)

integer,intent(in)::pk
real*8,intent(in)::t
real*8::ex

if(testType==0)then
 ex=0. 
else if(testType==1)then
 ex=exp(-((px(pk)-0.1)/(0.2*dp(pk)))**2.)
else if(testType==2)then
 ex=exp(-((px(pk)-0.1)/(0.2*dp(pk)))**2.)
else if(testType==3)then
 if(t<tol)then
  ex=exp(-((px(pk)-0.1)/(0.2*dp(pk)))**2.)
 else
  ex=0.
 end if
end if

!ex=300.*ex/dt

end function E


end program FokkerPlanck_CC
