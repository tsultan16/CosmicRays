module constants_mod
implicit none

!-----------------------------------------------------------------------------------------
integer,parameter::nt=100000!000
integer,parameter::nx=512
integer,parameter::ny=512
real*8,parameter::xmin=0.0,xmax=1.0,ymin=0.0,ymax=1.0,cour=0.8

integer,parameter::debug=0 !0:off 1:on
integer,parameter::discontinuityDir=1!1:x, 2:y, 3:oblique, 4:partial grid circle 
integer,parameter::boundaryType=4 !1:outflow 2:periodic in x, outflow in y, 3:reflecting 4:periodic, 5:peridoc in x reflecting in y
integer,parameter::initOption=0
integer,parameter::fluxType=3!1:Roe 2:HLL 3:HLLI
integer,parameter::HLLOption=0

integer,parameter::HLLIentropymode=1
integer,parameter::HLLIslowmode=1
integer,parameter::HLLIalfvenmode=1
integer,parameter::HLLIfastmode=0

integer,parameter::gravity=1 !0:off 1:on , uniform external gravitational field

integer,parameter::perturbationType=2!1:Sinusoidal 2:Random
integer,parameter::KH_test=0 !0:off 1:on
integer,parameter::RT_test=2 !0:off 1: vertical 2:radial
integer,parameter::cloudShock=0 !0:off 1:on
real,parameter::eps=0.00001

real*8,parameter::tol=1.d-20
real*8,parameter::min_pres=1.d-5
real*8,parameter::min_dens=1.d-5
real*8,parameter::pi=3.14159265359

real*8::gx,gy !gravitaional field components
real*8::divU(-1:nx+2,-1:ny+2)

!-------------------------------------------------------------------------------------------
integer,parameter::tracerType=2 !1:Mass Flux 2: Velocity Field
integer::trdebug=0 !0:off, 1:on
integer,parameter::init_option=2! 1:match fluid density distribution 2:uniformly distribute tracers in cells [minCell,maxCell] 3:Match fluid distrinution in cells [minCell,maxCell] 
integer,parameter::trboundaryType=2 !1:periodic, 2:outflow, 3:periodic in x & outflow in y
integer,parameter::N=1 !fixed # of tracers 
integer,parameter::minCellx=75,maxCellx=75,minCelly=75,maxCelly=75 
integer,parameter::interpolation_option=1 ;!nearest_cell , 2:(slope limited)linear interpolation
!-------------------------------------------------------------------------------------------
integer,parameter::tSkip=1000
!-------------------------------------------------------------------------------------------


end module constants_mod
