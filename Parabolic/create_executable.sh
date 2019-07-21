#!/bin/sh
gfortran ParabolicFokkerPlanck.f90 -o parabolic -ffpe-trap=invalid,zero,overflow,underflow
./parabolic
