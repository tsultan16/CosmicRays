module constants_mod
implicit none

!-----------------------------------------------------------------------------------------
integer,parameter::nt=1500
integer,parameter::nx=500
integer,parameter::np=100
integer,parameter::debug=0!0:off 1:on
integer,parameter::boundaryType=4 !1:outflow 2:periodic in x, outflow in y, 3:reflecting 4:periodic, 5:peridoc in x reflecting in y

real*8,parameter::xmin=0.0,xmax=5.0,cour=0.8
real*8,parameter::Q_user=2.0 
real*8::premin=1.d-25,densmin=1.d-25
real*8,parameter::tol=1.d-15

!----------------------------------------------------------------
real*8::u1(-2:nx+1,3),u2(-2:nx+1,3),flux(-2:nx+1,3)
real*8::u3(-1:nx-1,3)
real*8::ustar(-1:nx-1)
real*8::dt,dx !time step and spatial resolution
real*8::x
real*8::gam,gamil,gamul,gamel,gamee,gamuu
real*8::lambda(-1:nx-1,3)
!----------------------------------------------------------------
!----------------------------------------------------------------
integer,parameter::momentumMeshType=2 !1:uniform, 2:logarithmic
integer,parameter::testType=0 !0,1,2,3
integer,parameter::methodType=1 !1:C-C Fully Implicit, 2:C-C Semi-Implicit (Crank-Nicholson) 3:Simple Explicit, 4:Simple Implicit

real*8,parameter::cl=100!speed of light
real*8::Pc(-1:nx+2),dudx(0:nx-1)

real*8::f(-1:np,-1:nx),ftemp1(-1:np),ftemp2(-1:nx)
real*8::fp(0:np-1),fx(0:nx-1)
real*8::At(0:np-1),Bt(0:np-1),Ct(0:np-1)
real*8::Atx(0:nx-1),Btx(0:nx-1),Ctx(0:nx-1)
real*8::px(-2:np+3),dp(-1:np+2),dp2(-1:np+2)
real*8::pmin,pmax,t,nparticle

integer::shockCell(1)

integer::funit
!-----------------------------------------------------------------

integer,parameter::tSkip=10

end module constants_mod
