import numpy as np
import pandas as pd
import sys
import os
import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as unit

sys.path.append('../..')

from OpenSMOG3SPN2.forcefields.parsers import DNA3SPN2Parser
from OpenSMOG3SPN2.forcefields import SMOG3SPN2Model
from OpenSMOG3SPN2.utils.helper_functions import get_WC_paired_seq
from OpenSMOG3SPN2.utils.insert import insert_molecules

box_a, box_b, box_c = 20.0, 20.0, 20.0
n_dsDNA = 5

with open('single_nucl_dna_seq.txt', 'r') as f:
    seq = f.readlines()[0].strip()
full_seq = seq + get_WC_paired_seq(seq)

dna_parser = DNA3SPN2Parser.from_atomistic_pdb('dna.pdb', 'cg_dna.pdb', new_sequence=full_seq, temp_name='dna')
if not os.path.exists('start.pdb'):
    insert_molecules('cg_dna.pdb', 'start.pdb', n_mol=n_dsDNA, box=[box_a, box_b, box_c])

top = app.PDBFile('start.pdb').getTopology()
init_coord = app.PDBFile('start.pdb').getPositions()
model = SMOG3SPN2Model()
for i in range(n_dsDNA):
    model.append_mol(dna_parser)

model.create_system(top, box_a=box_a, box_b=box_b, box_c=box_c)
model.add_dna_bonds(force_group=5)
model.add_dna_angles(force_group=6)
model.add_dna_stackings(force_group=7)
model.add_dna_dihedrals(force_group=8)
model.add_dna_base_pairs(force_group=9)
model.add_dna_cross_stackings(force_group=10)
model.parse_all_exclusions()
model.add_all_vdwl(force_group=11)
model.add_all_elec(force_group=12)
temperature = 300*unit.kelvin
friction_coeff = 0.01/unit.picosecond
timestep = 10*unit.femtosecond
integrator = mm.LangevinMiddleIntegrator(temperature, friction_coeff, timestep)
model.set_simulation(integrator, platform_name='CUDA', init_coord=init_coord)
model.simulation.minimizeEnergy()
model.add_reporters(report_interval=10000, output_dcd='traj.dcd')
model.simulation.context.setVelocitiesToTemperature(temperature)
model.simulation.step(100000)



