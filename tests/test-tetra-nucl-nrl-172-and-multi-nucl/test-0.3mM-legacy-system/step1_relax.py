import numpy as np
import pandas as pd
import sys
import os
import argparse
import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as unit
import mdtraj
import time

sys.path.append('../../..')
from OpenSMOG3SPN2.forcefields.rigid import createRigidBodies
from OpenSMOG3SPN2.utils.chromatin_helper_functions import get_chromatin_rigid_bodies

parser = argparse.ArgumentParser()
parser.add_argument('--output_dcd', default='relax.dcd', help='output dcd path')
parser.add_argument('--output_interval', type=int, default=100000, help='output interval')
parser.add_argument('--n_steps', type=int, default=1000000, help='simulation number of steps')
args = parser.parse_args()

build_dir = '../build-0.3mM-legacy-system'
nonrigid_system_xml = f'{build_dir}/nonrigid_system.xml'
pdb = f'{build_dir}/start.pdb'

top = app.PDBFile(pdb).getTopology()
n_atoms = top.getNumAtoms()
n_atoms_per_single_nucl = app.PDBFile(f'{build_dir}/cg_single_nucl.pdb').getTopology().getNumAtoms()
init_coord = app.PDBFile(pdb).getPositions()
rigid_coord = init_coord

# set rigid bodies
# fix the whole tetranucleosome, and for each single nucleosome, fix histone core with middle 73 bp core DNA
rigid_bodies = []
n_tetra_nucl_atoms = app.PDBFile(f'{build_dir}/cg_tetra_nucl.pdb').getTopology().getNumAtoms()
#print(f'{n_tetra_nucl_atoms} atoms in tetranucleosome.')
rigid_bodies.append(list(range(n_tetra_nucl_atoms)))

single_nucl_rigid_body = np.array(get_chromatin_rigid_bodies(n_nucl=1, nrl=147)[0])
n_single_nucl = int((n_atoms - n_tetra_nucl_atoms)/n_atoms_per_single_nucl)
print(f'{n_single_nucl} single nucleosomes in the system.')
for i in range(n_single_nucl):
    rigid_bodies.append((single_nucl_rigid_body + n_tetra_nucl_atoms + i*n_atoms_per_single_nucl).tolist())

with open(nonrigid_system_xml, 'r') as f:
    system = mm.XmlSerializer.deserialize(f.read())
createRigidBodies(system, rigid_coord, rigid_bodies)

temperature = 300*unit.kelvin
friction_coeff = 0.01/unit.picosecond
timestep = 10*unit.femtosecond
integrator = mm.LangevinMiddleIntegrator(temperature, friction_coeff, timestep)
platform_name = 'CUDA'
platform = mm.Platform.getPlatformByName(platform_name)
properties = {'Precision': 'mixed'}
simulation = app.Simulation(top, system, integrator, platform, properties)
simulation.context.setPositions(init_coord)
simulation.minimizeEnergy()
simulation.context.setVelocitiesToTemperature(temperature)

dcd_reporter = app.DCDReporter(args.output_dcd, args.output_interval, enforcePeriodicBox=True)
state_reporter = app.StateDataReporter(sys.stdout, args.output_interval, step=True, time=True, potentialEnergy=True,
                                       kineticEnergy=True, totalEnergy=True, temperature=True, speed=True)
simulation.reporters.append(dcd_reporter)
simulation.reporters.append(state_reporter)
simulation.step(args.n_steps)


