#!/bin/sh
gfortran -c constants.f90
gfortran -c read_tracer.f90
gfortran -c CRtracer.f90
gfortran -c ParabolicFokkerPlanck.f90 
gfortran constants.o  read_tracer.o ParabolicFokkerPlanck.o CRtracer.f90 -o CRtracer -ffpe-trap=invalid,zero,overflow

