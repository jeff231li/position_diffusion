from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import netCDF4 as netcdf
import time
import sys
import numpy as np

dt       = 1.0
kr       = 30.0
Temp     = 300.0
Ntherm   = 50000
Nequil   = 500000
NLog     = 5000
Filename = 'Na-Box-30'

pos_ion  = np.zeros((Nequil,3), dtype=np.float64)
t_ion    = np.linspace(1,Nequil,Nequil)

# Hardware
sys.stdout = open('equilibrate.log', 'w')
platform   = Platform.getPlatformByName('CUDA')
properties = {'CudaDeviceIndex': '0', 'CudaPrecision': 'mixed'}

# Load CHARMM files
psf = CharmmPsfFile(Filename+'.psf')
pdb = PDBFile(Filename+'.pdb')
params = CharmmParameterSet('toppar_water_ions.str')
psf.setBox(3.0,3.0,3.0)

# Create System Object
system = psf.createSystem(params, 
                          nonbondedMethod=PME,
                          nonbondedCutoff=1.2*nanometer,
                          switchDistance=1.0*nanometer,
                          constraints=HBonds,
                          rigidWater=True,
                          removeCMMotion=True)

# Restrain Ion at the centre
center_of_mass = openmm.CustomExternalForce('(kr/2.0)*((x-x0)^2 + (y-y0)^2 + (z-z0)^2)')
center_of_mass.addGlobalParameter("kr", kr*kilocalories_per_mole/angstrom**2)
center_of_mass.addGlobalParameter("x0", 15*angstrom)
center_of_mass.addGlobalParameter("y0", 15*angstrom)
center_of_mass.addGlobalParameter("z0", 15*angstrom)
center_of_mass.addParticle(0, [])
system.addForce(center_of_mass)

## THERMOSTAT & BAROSTAT
#-------------------------------------------------------------------------------#
integrator = LangevinIntegrator(Temp*kelvin,      # Temperature of head bath
                                1.0/picosecond,   # Friction coefficient
                                dt*femtoseconds)  # Time step

barostat = MonteCarloBarostat(1.0*bar,Temp*kelvin,25)
system.addForce(barostat)

## INITIAL EQUILBRATION
#-------------------------------------------------------------------------------#
simulation = Simulation(psf.topology, system, integrator, platform, properties)
simulation.context.setPositions(pdb.positions)

print("Minimizing Energy")
sys.stdout.flush()
simulation.minimizeEnergy(maxIterations=1000, tolerance=5*kilojoule/mole)
simulation.context.setVelocitiesToTemperature(Temp*kelvin)
print(simulation.context.getState(getEnergy=True).getPotentialEnergy())

print("\nThermalizing for %i steps" % Ntherm)
sys.stdout.flush()
simulation.reporters.append(StateDataReporter('openmm-run.log',
                            NLog,
                            step=True,
                            kineticEnergy=True,potentialEnergy=True,totalEnergy=True,
                            temperature=True,
                            totalSteps=Ntherm+Nequil,
                            progress=True,remainingTime=True,speed=True,
                            separator="\t"))
simulation.step(Ntherm)
print(simulation.context.getState(getEnergy=True).getPotentialEnergy())

print("\nSaving Coordinates")
sys.stdout.flush()
positions = simulation.context.getState(getPositions=True).getPositions()
PDBFile.writeFile(simulation.topology, positions, open('thermalized.pdb', 'w'))

print("Running Production: %i" % Nequil)
sys.stdout.flush()
start = time.time()

for i in range(Nequil):
	simulation.step(1)
	pos_ion[i,:] = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)[0,:]*10

np.savetxt('Ion-COM-X-%ikcal.txt' % kr, np.c_[t_ion,pos_ion[:,0]], fmt="%10.5f")
np.savetxt('Ion-COM-Y-%ikcal.txt' % kr, np.c_[t_ion,pos_ion[:,1]], fmt="%10.5f")
np.savetxt('Ion-COM-Z-%ikcal.txt' % kr, np.c_[t_ion,pos_ion[:,2]], fmt="%10.5f")

end = time.time()
print("\tTotal time taken: %.2f sec" % (end-start))
