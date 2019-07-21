!Chang-Cooper(1970) algorithm for Isotropic Cosmic Ray-Fokker Planck Equation



program FokkerPlanck_Parabolic

implicit none

integer,parameter::np=100 !momentum mesh size
integer,parameter::nx=100 !spatial mesh size
integer,parameter::nt=200 !time steps
integer,parameter::meshType=2 !1:uniform, 2:logarithmic
integer,parameter::testType=0 !0,1,2,3
integer,parameter::methodType=1 !1:Fully Implicit, 2:Semi-Implicit (Crank-Nicholson) 3:Fully Explicit
integer,parameter::tSkip=1

real,parameter::tol=1.D-8
real*8::f(-1:np+2,-1:nx+2),ftemp1(-1:np),ftemp2(-1:nx)
real*8::fp(0:np-1),fx(0:nx-1)
real*8::At(0:np-1),Bt(0:np-1),Ct(0:np-1)
real*8::Atx(0:nx-1),Btx(0:nx-1),Ctx(0:nx-1)
real*8::px(-2:np+3),dp(-1:np+2),dp2(-1:np+2)
real*8::pmin,pmax,t,nparticle
real*8::xmin,xmax,dx
real*8::dt

integer::i,j,k,funit

if(methodType==1)then
 funit=10
 open(unit=funit,file='fp.txt')
 print*,'Cooper-Chang Fully Implicit Fokker-Planck Solver.'
end if

if(methodType==2)then
 funit=11
 open(unit=funit,file='fp2.txt')
 print*,'Cooper-Chang Semi Implicit Fokker-Planck Solver.'
end if

if(methodType==3)then
 funit=12
 open(unit=funit,file='fp3.txt')
 print*,'Cooper-Chang Explicit Fokker-Planck Solver.'
end if


open(unit=13,file='fx.txt')
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
do i=-2,np+3
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
 !print*,'i,p,dp=',i,px(i),dp(i)
end do

!Time Step Size
dt=1000.
do i=1,np-1
 dt=min(dt,(0.45)*dx**2./Kxx(px(i),0._8))
end do

dt=10.*dt

print*,'dt=',dt

do i=0,np-1
 do j=0,nx-1
  !Initial profiles
  f(i,j)=exp(-( (xmin+j*dx-0.5)/(10.*dx) )**2. )
 end do
 f(-1,j)=f(0,j)  
 f(np,j)=f(np-1,j)
end do


fp=0.
do i=0,np-1
  do j=0,nx-1
  fp(i)=fp(i)+dx*f(i,j)
 end do
  write(funit,*) px(i),fp(i)!+1.D-80
 !print*,'p,f(p)',px(i),fp(i)
end do

fx=0.
do j=0,nx-1
  do i=0,np-1
  fx(j)=fx(j)+dp(i)*f(i,j)
 end do
 write(13,*) xmin+j*dx,fx(j)!+1.D-50
 !print*,'x,f(x)',xmin+j*dx,fx(j)
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
  do j=0,nx-1
   ftemp1=0.
   do k=-1,np
    ftemp1(k)=f(k,j)
   end do
    
    if(methodType==3)then
      call computeFlux(ftemp1)  

    else
      call computeMomentumCoefficients()
      call tridiagSolve(np,ftemp1)
      
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
    
   call computeSpatialCoefficients(j)
   call tridiagSolveSpace(nx,ftemp2)
      
   ftemp2(-1)=ftemp2(0)  
   ftemp2(np)=ftemp2(np-1)

   do k=-1,nx
    f(j,k)=ftemp2(k)   
   end do
  end do



  fp=0.
  do k=0,np-1
   do j=0,nx-1
   fp(k)=fp(k)+dx*f(k,j)
   end do
   if(mod(i,tSkip)==0)then
    write(funit,*) px(k),fp(k)!+1.d-80
   end if
  end do

  fx=0.
  do j=0,nx-1
   do k=-1,np
   fx(j)=fx(j)+dp(k)*f(k,j)
   end do
   if(mod(i,tSkip)==0)then
    write(13,*) xmin+j*dx,fx(j)
   end if
  end do
  

  t=t+i*dt
  print*,'% Complete=',(real(i)/real(nt))*100.
  !if(i<0.1*nt)dt=1.2*dt

end do

close(unit=funit)
close(unit=13)

print*,'DONE.'

contains

subroutine computeFlux(ftemp1)
!Input Variables
real*8::ftemp1(-1:np)
!Local Variables
integer::jj,fupdate(-1:np),pflux(-1:np+2)
real*8::w(-2:np+2),Wplus(-1:np+1),Wminus(-1:np+1),Wbig(-1:np+1)
real*8::B2,C2,B1,C1,pi2,pi2L

!-------------------------------------------------------------
!Weights: w_j, W_j, W^+_j+1/2, W^-_j+1/2 , F^n_i+1/2
!-------------------------------------------------------------
do jj=-1,np+2
 B2=B(px(jj),0._8,0._8)
 C2=C(px(jj),0._8)
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
 C2=C(px(jj),0._8)
 pflux(jj)=(C2/dp2(jj))*(WPlus(jj)*ftemp1(jj+1)-WMinus(jj)*ftemp1(jj))
end do

!No-flux Boundary Condition
pflux(-1)=0.
pflux(np-1)=0.

!Update distribution function
!fupdate=ftemp1
do jj=0,np-1
ftemp1(jj)=ftemp1(jj)*(1.+dt*Hp(px(jj),0._8)) &
           +(dt/A(px(jj)))*(pflux(jj)-pflux(jj-1))/dp(jj) &
           + dt*S(jj,t)
end do
!ftemp1=fupdate

end subroutine computeFlux

subroutine computeMomentumCoefficients()

!Input variables

integer::jj
real*8::w(-2:np+2),Wplus(-1:np+1),Wminus(-1:np+1),Wbig(-1:np+1)
real*8::B2,C2,B1,C1,pi2,pi2L

!-------------------------------------------------------------
!Weights: w_j, W_j, W^+_j+1/2, W^-_j+1/2 
!-------------------------------------------------------------
do jj=-1,np+2
 B2=B(px(jj),0._8,0._8)
 C2=C(px(jj),0._8)


 w(jj)=dp(jj)*B2/C2 

 if(abs(w(jj))>1.d4 .or. abs(w(jj))<1.d-12)then
   print*,'***ERROR***, w_j=',w(jj)
   print*,'dp,B2,C2=',dp(jj),B2,C2
   print*,'w_j is out of bounds.' 
   print*,'Need to re adjust size of momentum advection and diffusion co-efficients.' 
   STOP
 end if


 !print*,'j,dp,w=',j,w(jj)  
 !if(abs(w(jj))<0.1)then
  ! Wbig(jj)=1./(1.+(1./24.)*w(jj)**2.+(1./1920.)*w(jj)**4.)
 !else
  !Wbig(jj)=abs(w(jj))*exp(-0.5*abs(w(jj)))/(1.-exp(-abs(w(jj))))
 !end if

  Wbig(jj)=0.5*w(jj)/sinh(0.5*w(jj))
  !print*,'j,w,Wbig=',j,w(j),Wbig(j)
end do

do jj=-1,np+1
 !print*,'Wbig=',Wbig(j)
 !if(w(jj)<tol) then
 !  Wplus(jj)=0.5*(Wbig(jj)*exp(0.5*w(jj))+Wbig(jj+1)*exp(0.5*w(jj+1)))
 !  Wminus(jj)=0.5*(Wplus(jj)-w(jj)+Wplus(jj+1)-w(jj+1))  
 !else
   Wminus(jj)=0.5*(Wbig(jj)*exp(-0.5*w(jj))+Wbig(jj+1)*exp(-0.5*w(jj+1)) )
   Wplus(jj)=0.5*(Wbig(jj)*exp(0.5*w(jj))+Wbig(jj+1)*exp(0.5*w(jj+1)) )
 !end if
 !Wplus(jj)=1./(1.-exp(-w(jj)))
 !Wplus(jj)=Wplus(jj)+1./(1.-exp(-w(jj+1)))
 !Wplus(jj)=0.5*Wplus(jj)
 !Wminus(jj)=1./(exp(w(jj))-1.)
 !Wminus(jj)=Wminus(jj)+1./(exp(w(jj+1))-1.)
 !Wminus(jj)=-0.5*Wminus(jj)
 !print*,'jj,Wplus,Wminus=',jj,Wplus(jj),Wminus(jj)
end do


!No-flux Boundary condition requirement
Wplus(-1)=0.
Wminus(-1)=0.
Wplus(np-1)=0.
Wminus(np-1)=0.

!At=0.
!Bt=0.
!Ct=0.

!------------------------------------------------------------
!Co-efficients for tridiagonal system: ~A_j,~B_j,~C_j
!------------------------------------------------------------
do jj=0,np-1
  C1=0.5*( C(px(jj-1),0._8)+C(px(jj),0._8) )  !C_j-1/2
  C2=0.5*( C(px(jj),0._8)+C(px(jj+1),0._8) )  !C_j+1/2

  At(jj)=dt/(dp(jj)*dp2(jj-1)*A(px(jj)))
  At(jj)=At(jj)*C1*Wminus(jj-1)
  At(jj)=-At(jj)
  if(methodType==2)then
   At(jj)=At(jj)*0.5
  end if

  Bt(jj)= C1*Wplus(jj-1)/dp2(jj-1) &
        +C2*Wminus(jj)/dp2(jj)
  Bt(jj)=Bt(jj)*dt/(dp(jj)*A(px(jj)))
  Bt(jj)=Bt(jj)+1.+dt*Hp(px(jj),0._8)
  if(methodType==2)then 
   Bt(jj)=0.5*(Bt(jj)+1.)
  end if

  Ct(jj)=dt/(dp(jj)*dp2(jj)*A(px(jj)))
  Ct(jj)=Ct(jj)*C2*Wplus(jj)
  Ct(jj)=-Ct(jj)
  if(methodType==2)then
   Ct(jj)=Ct(jj)*0.5
  end if

  !print*,'jj,At,Bt,Ct=',jj,At(jj),Bt(jj),Ct(jj)
end do
 
end subroutine computeMomentumCoefficients


subroutine computeSpatialCoefficients(pk)

!Input variables
integer::pk
integer::jj
real*8::w(-2:nx+2),Wplus(-1:nx+1),Wminus(-1:nx+1),Wbig(-1:nx+1)
real*8::B2,C2,B1,C1,pi2,pi2L

!-------------------------------------------------------------
!Weights: w_j, W_j, W^+_j+1/2, W^-_j+1/2 
!-------------------------------------------------------------
do jj=-1,nx+2
 B2=E(px(pk),0._8,0._8)
 C2=Kxx(px(pk),0._8)
 w(jj)=dx*B2/C2 

 if(abs(w(jj))>1.d4 .or. abs(w(jj))<1.d-12)then
   print*,'***ERROR***, w_j=',w(jj)
   print*,'dp,B2,C2=',dp(pk),B2,C2
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


!------------------------------------------------------------
!Co-efficients for tridiagonal system: ~A_j,~B_j,~C_j
!------------------------------------------------------------
do jj=0,nx-1
  C1=Kxx(px(pk),0._8)
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
  Btx(jj)=Btx(jj)+1.+dt*Hx(px(pk),0._8,0._8)
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


subroutine tridiagSolve(n,ftemp) !from Numerical Recipes book

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

end subroutine tridiagSolve

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
function A(p) result(ax)
real*8,intent(in)::p
real*8::ax
  
 ax=1.!p*p

end function A

!Momentum Advection co-efficient
function B(p,t,dudx) result(bx)
real*8,intent(in)::p,t,dudx
real*8::bx

 bx=0.000001!1.!(p*p)!(p**3.)*dudx/3.

end function B

!Momentum Diffusion co-efficient
function C(p,t) result(cx)

real*8,intent(in)::p,t
real*8::cx,Dpp

 Dpp=0.00001
 cx=Dpp*p*p

end function C

!Spatial Advection co-efficient
function E(p,t,dudx) result(ex)
real*8,intent(in)::p,t,dudx
real*8::ex

ex=0.1

end function E

!Spatial Diffusion co-efficient
function Kxx(p,t) result(kx)

real*8,intent(in)::p,t
real*8::kx

kx=0.01!*(p**0.6)

end function Kxx


!Momentum Source term co-efficient
function Hp(dudx,t) result(hx)
 
real*8,intent(in)::dudx,t
real*8::hx

hx=0.!-dudx/3.

end function Hp


!Spatial Source term co-efficient
function Hx(p,dudx,t) result(hxx)
 
real*8,intent(in)::p,dudx,t
real*8::hxx

hxx=dudx

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

lx=0.

end function L




subroutine thomasAlgorithm(n,ftemp)
!Input variables
integer::n
real*8::ftemp(-1:n)

integer::kk
real*8::e(0:n-1),g(0:n-1)


!------------------------------------------------------------
!Compute Thomas algorithm co-efficients: e_k,g_k
!------------------------------------------------------------ 
e(0)=0.
g(0)=0. !These correspond to the boundary condition f(0)=0.

do kk=0,n-1
 e(kk)=Ct(kk)/(Bt(kk)-e(kk-1)*At(kk))
 g(kk)=(ftemp(kk)+At(kk)*g(kk-1))/(Bt(kk)-e(kk-1)*At(kk))
end do

do kk=n-1,0,-1
  ftemp(kk)=g(kk)+e(kk)*ftemp(kk+1)
end do

end subroutine thomasAlgorithm


end program FokkerPlanck_Parabolic
