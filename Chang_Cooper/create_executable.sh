#!/bin/sh
gfortran changCooper.f90 -o CC -ffpe-trap=invalid,zero,overflow
./CC
