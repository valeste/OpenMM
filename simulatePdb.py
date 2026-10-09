
from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
from pathlib import Path

# Define output folder within your OpenMM project
output_dir = Path.home() / "Documents" / "OpenMM" / "output"
output_dir.mkdir(parents=True, exist_ok=True)

pdb = PDBFile('input.pdb')
forcefield = ForceField('amber19-all.xml', 'amber19/tip3pfb.xml')
system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)
simulation.minimizeEnergy()
# Below we add a reporter to the simulation, which will write out the positions of the atoms in the system to a DCD file every 1000 steps.
simulation.reporters.append(DCDReporter(str(output_dir / 'output.dcd'), 1000))
# Below we add a reporter to the simulation, which will print out information about the simulation in the terminal window to the console every 1000 steps.
# The information includes the current step number, the potential energy of the system, and the temperature of the system.
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))
# Save temperature and potential energy as a CSV
simulation.reporters.append(
    StateDataReporter(
        str(output_dir / 'simulation_data.csv'),
        1000,
        step=True,
        potentialEnergy=True,
        temperature=True
    )
)
simulation.step(10000)