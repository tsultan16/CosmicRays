!Chang-Cooper(1970) algorithm for Cosmic Ray-Fokker Planck Equation
!(evolution in momentum space only+isotropic plasma)


program FokkerPlanck_Parabolic

implicit none

integer,parameter::np=100 !momentum mesh size
integer,parameter::nx=10 !spatial mesh size
integer,parameter::nt=600 !time steps
integer,parameter::meshType=1 !1:uniform, 2:logarithmic
integer,parameter::testType=0 !0,1,2,3
integer,parameter::methodType=1 !1:Fully Implicit, 2:Semi-Implicit (Crank-Nicholson)

real,parameter::tol=1.D-8
real*8::f(-1:np+2,-1:nx+2),ftemp1(0:np),ftemp2(0:nx)
real*8::At(0:np),Bt(0:np),Ct(0:np)
real*8::px(-2:np+3),dp(0:np+1),dp2(0:np+1)
real*8::pmin,pmax,t,nparticle
real*8::xmin,xmax,dx
real*8::dt=100.

integer::i,j,k

if(methodType==1)then
 open(unit=10,file='fp.txt')
 print*,'Cooper-Chang Fully Implicit Fokker-Planck Solver.'
else if(methodType==2)then
 open(unit=10,file='fp2.txt')
 print*,'Cooper-Chang Semi Implicit Fokker-Planck Solver.'
end if
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

nparticle=0.

!Generate momentum mesh points
do i=-1,np+2
  if(meshType==1)then
   px(i)=pmin+i*(pmax-pmin)/np
  else if(meshType==2)then
   px(i)=pmin*(pmax/pmin)**(real(i)/real(np))
  end if
end do

!Momentum mesh spacing
do i=-1,np+2
 dp(i)=0.5*(px(i+1)-px(i-1))
 dp2(i)=px(i+1)-px(i)
end do

if(testType==3)then
 do i=0,np
  nparticle=nparticle+exp(-((px(i)-0.1)/(0.1*dp(i)))**2.)*dp(i)
 end do
else
 do i=0,np
  nparticle=nparticle+exp(-((px(i)-1.)/(0.05*dp(i)))**2.)*dp(i)
 end do
end if
 
do i=0,np
  do j=-1,nx+1
  !Initial profile
  !if(testType==3)then 
    f(i,j)=exp(-((px(i)-500.)/(0.5*dp(i)))**2.)/nparticle 
  !else if(testType==0)then
  !  f(i,j)=1.D-15*px(i)**(-5.) 
  !else
    !f(i,j)=0.+1.D-60 
  !end if

end do

f(-1,j)=0.
f(np+1,j)=0.
end do


do i=1,np
  write(10,*) px(i),f(i,1)+1.D-50
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

   do k=0,np
    ftemp1(k)=f(k,j)
   end do

   call computeMomentumCoefficients(j)
   call tridiagSolve(np,ftemp1)
   
   do k=0,np
    f(k,j)=ftemp1(k)
   end do

  end do

  !-------------------------------------
  !Spatial Operator Sweeps   
  !------------------------------------- 
  do j=0,np

  end do


  do j=1,np
   !print*,'p,f(p)=',px(j),u(j)
    write(10,*) px(j),f(j,1)+1.D-50 
  end do


  t=t+i*dt
  print*,'% Complete=',(real(i)/real(nt))*100.


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
!Weights: w_j+1/2, W_j+1/2, W^+_j+1/2, W^-_j+1/2 
!-------------------------------------------------------------
do j=-1,np+2
 B2=0.5*( B(px(j),0._8,0._8)+B(px(j+1),0._8,0._8) )
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
  if(methodType==2)then
   At(j)=At(j)*0.5
  end if

  Bt(j)= C1*Wplus(j-1)/dp2(j-1) &
        +C2*Wminus(j)/dp2(j)
  Bt(j)=Bt(j)*dt/(dp(j)*A(px(j)))
  Bt(j)=Bt(j)+1.-dt*Hp(px(j),0._8)
  if(methodType==2)then 
   Bt(j)=0.5*(Bt(j)+1.)
  end if

  Ct(j)=dt/(dp(j)*dp2(j)*A(px(j)))
  Ct(j)=Ct(j)*C2*Wplus(j)
  Ct(j)=-Ct(j)
  if(methodType==2)then
   Ct(j)=Ct(j)*0.5
  end if
end do

!Impose No-flux boundary conditions
At(0)=0.
Ct(np)=0. 

end subroutine computeMomentumCoefficients

subroutine tridiagSolve(n,ftemp) !from Numerical Recipes book

!Input variables
integer::n
real*8::ftemp(0:n)

real*8::bet,gam(0:n),r(0:n)
integer::k

do k=0,np
 r(k)=ftemp(k)+dt*S(k,t)
 if(methodType==2)then
   r(k)=r(k)+ftemp(k)-At(k)*ftemp(k-1)-Bt(k)*ftemp(k)-Ct(k)*ftemp(k+1)
 end if
end do

bet=Bt(0)
ftemp(0)=r(0)/bet

!Upper triangular decomposition and forward substitution
do k=1,n
  gam(k)=Ct(k-1)/bet
  bet=Bt(k)-At(k)*gam(k)
  ftemp(k)=(r(k)-At(k)*ftemp(k-1))/bet
end do

!Back substitution
do k=n-1,0,-1
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

 bx=0.D-15*(1.+p)!(p**3.)*dudx/3.

end function B

!Momentum Diffusion co-efficient
function C(p,t) result(cx)

real*8,intent(in)::p,t
real*8::cx,Dpp

 Dpp=1.D0
 cx=Dpp!*p*p

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
function Hp(dudx,t) result(hx)
 
real*8,intent(in)::dudx,t
real*8::hx

hx=0.!-1.!-dudx/3.

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

sx=0.!exp(-((px(pk)-0.1)/(0.05*dp(pk)))**2.)

end function S


end program FokkerPlanck_Parabolic
