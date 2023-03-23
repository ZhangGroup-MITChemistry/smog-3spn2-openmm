import numpy as np
import pandas as pd
import sys
import os
import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as unit

sys.path.append('../../..')

from OpenSMOG3SPN2.forcefields.parsers import SMOGParser
from OpenSMOG3SPN2.forcefields import SMOG3SPN2Model

histone_parser = SMOGParser.from_atomistic_pdb('../pdb-files/histone.pdb', 'histone_CA.pdb')

# parse dna later

nucl = SMOG3SPN2Model()
nucl.append_mol(histone_parser)
top = app.PDBFile('histone_CA.pdb').getTopology()
init_coord = app.PDBFile('histone_CA.pdb').getPositions()
nucl.create_system(top)
nucl.add_protein_bonds(force_group=1)
nucl.add_protein_angles(force_group=2)
nucl.add_protein_dihedrals(force_group=3)
nucl.add_native_pairs(force_group=4)

# add dna forces later

temperature = 300*unit.kelvin
friction_coeff = 0.01/unit.picosecond
timestep = 10*unit.femtosecond
integrator = mm.LangevinMiddleIntegrator(temperature, friction_coeff, timestep)
nucl.set_simulation(integrator, platform_name='CPU', init_coord=init_coord)
nucl.simulation.minimizeEnergy()
output_interval = 100
output_dcd = 'output.dcd'
nucl.add_reporters(output_interval, output_dcd)
nucl.simulation.context.setVelocitiesToTemperature(temperature)
nucl.simulation.step(500)



    
