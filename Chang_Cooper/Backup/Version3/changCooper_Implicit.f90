!Chang-Cooper(1970) algorithm for Cosmic Ray-Fokker Planck Equation
!(evolution in momentum space only+isotropic plasma)
!Fully-Implicit

program FokkerPlanck_CC

implicit none

integer,parameter::np=100 !momentum bins
integer,parameter::nt=500 !time steps
integer,parameter::meshType=2 !1:uniform, 2:logarithmic
integer,parameter::testType=3 !0,1,2,3
real,parameter::tol=1.D-8
real*8::u(-1:np+2),r(0:np)
real*8::delta(-1:np),At(0:np),Bt(0:np),Ct(0:np)
real*8::px(-2:np+3),dp(0:np+1),dp2(0:np+1)
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
do i=-1,np+2
 dp(i)=0.5*(px(i+1)-px(i-1))
 dp2(i)=px(i+1)-px(i)
end do


do i=0,np
  !Initial profile
  if(testType==3)then
    u(i)=exp(-((px(i)-0.1)/(0.1*dp(i)))**2.)*1.D5
  else
    u(i)=0.
  end if

  write(10,*) px(i),u(i)+1.D-50
     !added small number so that u not equal to zero, which is bad for log plots

end do

!----------------------------------------
!Evolve Distribution Function
!----------------------------------------
t=0._8
do i=1,nt
  print*,'Time Step=',i

  call computeCoefficients()
  
  call tridiagSolve()

  do j=1,np
   !print*,'p,f(p)=',px(j),u(j)
   write(10,*) px(j),u(j)+1.D-50 
      !added small number so that u not equal to zero, which is bad for log plots
  end do

  t=t+i*dt
  print*,'% Complete=',(real(i)/real(nt))*100.

end do

close(unit=10)

contains


subroutine computeCoefficients()
integer::j
real*8::w(-1:np+2),Wplus(-1:np+2),Wminus(-1:np+2),Wbig(-1:np+2)
real*8::B2,C2,B1,C1,pi2,pi2L

!-------------------------------------------------------------
!Weights: w_j+1/2, W_j+1/2, W^+_j+1/2, W^-_j+1/2 
!-------------------------------------------------------------
do j=-1,np+2
 B2=0.5*( B(px(j),0._8)+B(px(j+1),0._8) )
 C2=0.5*( C(px(j),0._8)+C(px(j+1),0._8) )
 !B2=B(0.5*(px(j)+px(j+1)),0._8)
 !C2=C(0.5*(px(j)+px(j+1)),0._8)
 !print*,'j,dp2,B2,C2=',j,dp2(j),B2,C2
 w(j)=dp2(j)*B2/C2
 !print*,'w=',w(j)
 if(abs(w(j))<0.1)then
   Wbig(j)=1./(1.+(1./24.)*w(j)**2.+(1./1920.)*w(j)**4.)
 else
   Wbig(j)=abs(w(j))*exp(-0.5*abs(w(j)))/(1.-exp(-abs(w(j))))
 end if
 !print*,'Wbig=',Wbig(j)
 if(w(j)<tol) then
   Wplus(j)=Wbig(j)*exp(0.5*w(j))
   Wminus(j)=Wplus(j)-w(j)  
 else
   Wminus(j)=Wbig(j)*exp(-0.5*w(j))
   Wplus(j)=Wminus(j)+w(j)
 end if 
 !print*,'Wplus,Wminus=',Wplus(j),Wminus(j)
end do


!------------------------------------------------------------
!Co-efficients for tridiagonal system: ~A_j,~B_j,~C_j
!------------------------------------------------------------
do j=0,np
  C1=0.5*( C(px(j-1),0._8)+C(px(j),0._8) )  !C_j+1/2
  C2=0.5*( C(px(j),0._8)+C(px(j+1),0._8) )  !C_j-1/2

  At(j)=dt/(dp(j)*dp2(j-1)*A(px(j)))
  At(j)=At(j)*C1*Wminus(j-1)
  At(j)=-At(j)

  Bt(j)= C1*Wplus(j-1)/dp2(j-1) &
        +C2*Wminus(j)/dp2(j)
  Bt(j)=Bt(j)*dt/(dp(j)*A(px(j)))
  Bt(j)=Bt(j)+1.-dt*D(px(j),0._8)

  Ct(j)=dt/(dp(j)*dp2(j)*A(px(j)))
  Ct(j)=Ct(j)*C2*Wplus(j)
  Ct(j)=-Ct(j)
end do

!Impose No-flux boundary conditions
At(0)=0.
Ct(np)=0. 

end subroutine computeCoefficients

subroutine tridiagSolve() !from Numerical Recipes book

real*8::bet,gam(0:np)
integer::k

do k=0,np
 r(k)=u(k)+dt*E(k,t)
end do

bet=Bt(0)
u(0)=r(0)/bet

!Upper triangular decomposition and forward substitution
do k=1,np
  gam(k)=Ct(k-1)/bet
  bet=Bt(k)-At(k)*gam(k)
  u(k)=(r(k)-At(k)*u(k-1))/bet
end do

!Back substitution
do k=np-1,0,-1
 u(k)=u(k)-gam(k+1)*u(k+1)
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
 ex=exp(-((px(pk)-0.1)/(0.05*dp(pk)))**2.)
else if(testType==2)then
 ex=exp(-((px(pk)-0.1)/(0.05*dp(pk)))**2.)
else if(testType==3)then
 ex=0.
end if

!ex=ex*1.D30

end function E


end program FokkerPlanck_CC
