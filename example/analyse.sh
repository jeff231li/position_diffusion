#!/bin/bash

gfortran -o ../diffusion ../diffusion.f90
../diffusion 2.0 0.001 Ion-COM-X-20kcal.txt ACF-X.txt
