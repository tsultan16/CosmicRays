!Chang-Cooper(1970) algorithm for Cosmic Ray-Fokker Planck Equation
!(evolution in momentum space only+isotropic plasma)


program FokkerPlanck_CC

implicit none

integer,parameter::np=200 !momentum bins
integer,parameter::nt=500 !time steps
real,parameter::tol=1.D-8
real*8::f(0:np+1),r(0:np+1)
real*8::delta(-1:np),At(0:np),Bt(0:np),Ct(0:np)
real*8::e(0:np-1),g(0:np-1)
real*8::pmin,pmax,dp,p
real*8::dt=100.

integer::i,j

open(unit=10,file='fp.txt')


!----------------------------------------------
!Initialization
!----------------------------------------------
!Momentum space bounds
pmin=0.0
pmax=1000.
dp=(pmax-pmin)/np


do i=1,np
  !Gaussian Initial Profile
  !f(i)=exp(-(((i-np/2)*dp)/(0.02*np*dp))**2.)
  !Square Initial Profile
  !if(i>0.4*np .and. i< 0.6*np)then
  !  f(i)=1.
  !else
  !  f(i)=0.
  !end if

  !Powerlaw Initial profile
   if(i>0.2*np .and. i<0.8*np)then
    f(i)=(i*dp)**(-5.)
   else
    f(i)=0.
  end if


  write(10,*) pmin+i*dp,f(i)

end do

!Boundary Conditions
f(1)=0.!f(1)
f(np)=0.!f(np)

!Evolve Distribution Function
do i=1,nt

  print*,'Time Step=',i
  call computeCoefficients()

  !f(np)=0.
  !do j=np-1,2,-1
  ! f(j)=f(j+1)*e(j)+g(j)
  !end do 

  r=f
  call tridiagSolve()

  f(0)=0.!f(1)
  f(np+1)=0.!f(np)
  
  do j=1,np
   write(10,*) pmin+j*dp,f(j)
  end do

  print*,'% Complete=',(real(i)/real(nt))*100.

end do

close(unit=10)

contains


subroutine computeCoefficients()
integer::k
real*8::w(np),Wj(np)

!--------------------------------------
!First compute weights: w_j, W_j
!--------------------------------------
do k=0,np
 w(k)=dp*B((k+0.5)*dp,0._8)/C((k+0.5)*dp,0._8)
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


!------------------------------------------------------------
!Compute co-efficients for tridiagonal system: ~A_j,~B_j,~C_j
!------------------------------------------------------------
do k=1,np
!wj=dp*B((k+0.5)*dp,0._8)/C((k+0.5)*dp,0._8)
!At(k)=dt*C((k+0.5)*dp,0._8)/(dp*dp*A(k*dp))
!At(k)=At(k)*((1.-delta(k))*wj+1.)
!At(k)=-At(k)

At(k)=dt/(dp*dp*A(k*dp))
At(k)=At(k)*C((k+0.5)*dp,0._8)*Wj(k)*exp(w(k))
At(k)=-At(k)


!Bt(k)=dt/(dp*A(k*dp))
!Bt(k)=Bt(k)*( (C((k+0.5)*dp,0._8)+C((k-0.5)*dp,0._8))/dp &
  !    +(1.-delta(k-1))*B((k-0.5)*dp,0._8)-delta(k)*B((k+0.5)*dp,0._8) )+1.

Bt(k)= C((k+0.5)*dp,0._8)*Wj(k) &
      +C((k-0.5)*dp,0._8)*Wj(k-1)*exp(w(k-1)) 
Bt(k)=Bt(k)*dt/(dp*dp*A(k*dp))+1.


!Ct(k)=dt/(dp*A(k*dp))
!Ct(k)=Ct(k)*( C((k-0.5)*dp,0._8)/dp-delta(k-1)*B((k-0.5)*dp,0._8)   )
!Ct(k)=-Ct(k)

Ct(k)=dt/(dp*dp*A(k*dp))
Ct(k)=Ct(k)*C((k-0.5)*dp,0._8)*Wj(k-1)
Ct(k)=-Ct(k)

end do

 
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


function A(p) result(ax)
real*8,intent(in)::p
real*8::ax
  
ax=0.01*p*p

end function A


function B(p,t) result(bx)
real*8,intent(in)::p,t
real*8::bx

bx=0.001 !Only diffusion for this test, no advection in momentum space

end function B


function C(p,t) result(cx)

real*8,intent(in)::p,t
real*8::cx

cx=1.*p !*p*p

end function C


end program FokkerPlanck_CC
