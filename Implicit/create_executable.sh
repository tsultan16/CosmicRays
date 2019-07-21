#!/bin/sh
gfortran ImplicitFokkerPlanck.f90 -o implicit -ffpe-trap=invalid,zero,overflow
./implicit
