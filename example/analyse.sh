#!/bin/bash

gfortran -O3 -o diffusion ../diffusion.f90
./diffusion 2.0 0.001 Ion-COMX.txt ACF-X.txt
