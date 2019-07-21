#!/bin/sh
gfortran -c constants.f90
gfortran -c RoeSolver_entropyfixed.f90
gfortran -c CGMVFokkerPlanck.f90 
gfortran constants.o RoeSolver_entropyfixed.o CGMVFokkerPlanck.o EulerMH.f90 -o EulerCR -ffpe-trap=invalid,zero,overflow
./EulerCR
