# position_diffusion
Position dependent diffusivity of a molecule under an external harmonic potential. The method is based on the derivation given in the reference below.

Hummer, Gerhard. "Position-dependent diffusion coefficients and free energies from Bayesian analysis of equilibrium and replica molecular dynamics simulations." New Journal of Physics 7.1 (2005): 34.

Compilation:
```
gfortran -o diffusion diffusion.f90
```
The code takes in four command line arguments:
```
1. Autocorrelation cutoff time (in ps)
2. Time step between frames (in ps)
3. Trajectory file name (two column ASCII)
4. Output name for ACF results
```
The trajectory should be formatted as a two column file with the time/frame number in the first and the position in the second. The positions should be in units of Angstroms. Running the program with
```
./diffusion 2.0 0.001 wind01.traj acf01.txt
```
should give an output file acf01.txt. The output file contains the some information in the header followed by the autocorrelation function.
```
# Frame time (ps)          :  0.00100
# Correlation time (ps)    :  2.00000
# Characteristic time (ps) :  0.24332
# Diffusion (x10^-9 m^2/s) :  1.25022
   0.00100   1.18259
   0.00200   1.18151
   0.00300   1.18012
   0.00400   1.17821
   0.00500   1.17590
   0.00600   1.17310
   0.00700   1.16977
   0.00800   1.16591
   ...
```
