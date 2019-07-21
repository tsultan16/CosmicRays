module constants_mod
implicit none

!-----------------------------------------------------------------------------------------
integer,parameter::nt=1000
integer,parameter::nx=100
integer,parameter::N=1 !# of tracers
integer,parameter::np=500 !# of momentum bins
integer,parameter::debug=0!0:off 1:on
integer,parameter::boundaryType=1 !1:outflow 2:periodic in x, outflow in y, 3:reflecting 4:periodic, 5:peridoc in x reflecting in y

real*8,parameter::xmin=0.0,xmax=1.0,ymin=0.0,ymax=1.0
real*8,parameter::tol=1.d-15
!-----------------------------------------------------------------------------------------
real*8::dt,dx !time step and spatial resolution
real*8::x,xd
!----------------------------------------------------------------
integer,parameter::momentumMeshType=2 !1:uniform, 2:logarithmic
integer,parameter::testType=0 !0,1,2,3
integer,parameter::methodType=1 !1:C-C Fully Implicit, 2:C-C Semi-Implicit (Crank-Nicholson) 

real*8,parameter::cl=100!speed of light

real*8::f(N,-1:np)
real*8::At(0:np-1),Bt(0:np-1),Ct(0:np-1)
real*8::Atx(0:nx-1),Btx(0:nx-1),Ctx(0:nx-1)
real*8::px(-2:np+3),dp(-1:np+2),dp2(-1:np+2)
real*8::pmin,pmax,t,nparticle

integer::shockCell(1)

integer::funit
!-----------------------------------------------------------------
!-----------------------------------------------------------------

integer,parameter::tSkip=10

end module constants_mod
