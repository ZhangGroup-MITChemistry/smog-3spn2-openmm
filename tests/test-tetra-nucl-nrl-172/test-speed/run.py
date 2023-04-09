import numpy as np
import pandas as pd
import sys
import os
import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as unit
import mdtraj

sys.path.append('../../..')

from OpenSMOG3SPN2.forcefields.parsers import SMOGParser, DNA3SPN2Parser
from OpenSMOG3SPN2.forcefields import SMOG3SPN2Model
from OpenSMOG3SPN2.utils.helper_functions import get_WC_paired_seq
from OpenSMOG3SPN2.forcefields.rigid import createRigidBodies
from OpenSMOG3SPN2.utils.chromatin_helper_functions import get_chromatin_rigid_bodies

n_nucl = 4
nrl = 172
platform_name = 'CUDA'
tetra_nucl = SMOG3SPN2Model()

# load proteins
for i in range(n_nucl):
    histone_i_parser = SMOGParser.from_atomistic_pdb(f'../pdb-files/histone_{i}.pdb', f'histone_{i}_CA.pdb', 
                                                     default_parse=False)
    histone_i_parser.parse_mol(get_native_pairs=False)
    tetra_nucl.append_mol(histone_i_parser)

# load DNA with the target sequence
# the sequence saved in dna_seq.txt only involves the first ssDNA
# we have to get the full WC-paired sequence
with open('dna_seq.txt', 'r') as f:
    seq1 = f.readlines()[0].strip()
seq2 = get_WC_paired_seq(seq1)
target_seq = seq1 + seq2

dna_parser = DNA3SPN2Parser.from_atomistic_pdb('../pdb-files/dna.pdb', 'cg_dna.pdb', new_sequence=target_seq)
tetra_nucl.append_mol(dna_parser)
tetra_nucl.atoms_to_pdb('cg_chromatin.pdb')
tetra_nucl.parse_all_exclusions()

top = app.PDBFile('cg_chromatin.pdb').getTopology()

# load coordinates
dcd = '../lammps-rerun/traj.dcd'
traj = mdtraj.load_dcd(dcd, top='cg_chromatin.pdb')
init_coord = traj.xyz[-1]*unit.nanometer
rigid_coord = init_coord
tetra_nucl.create_system(top, box_a=200, box_b=200, box_c=200)
rigid_bodies = get_chromatin_rigid_bodies(n_nucl, nrl)
createRigidBodies(tetra_nucl.system, rigid_coord, rigid_bodies)
#keep_unchanged_list = ['protein_exclusions', 'dna_exclusions', 'exclusions']
#tetra_nucl.set_rigid_bodies(rigid_coord, rigid_bodies, keep_unchanged=keep_unchanged_list)
tetra_nucl.add_protein_bonds(force_group=1)
tetra_nucl.add_protein_angles(force_group=2)
#tetra_nucl.add_protein_dihedrals(force_group=3)
#tetra_nucl.add_native_pairs(force_group=4)
tetra_nucl.add_dna_bonds(force_group=5)
tetra_nucl.add_dna_angles(force_group=6)
tetra_nucl.add_dna_stackings(force_group=7)
tetra_nucl.add_dna_dihedrals(force_group=8)
tetra_nucl.add_dna_base_pairs(force_group=9)
tetra_nucl.add_dna_cross_stackings(force_group=10)
tetra_nucl.add_all_vdwl(force_group=11)
tetra_nucl.add_all_elec(force_group=12)
#tetra_nucl.save_system('system.xml')

temperature = 300*unit.kelvin
friction_coeff = 0.01/unit.picosecond
timestep = 10*unit.femtosecond
integrator = mm.LangevinMiddleIntegrator(temperature, friction_coeff, timestep)
tetra_nucl.set_simulation(integrator, platform_name=platform_name, init_coord=init_coord)
tetra_nucl.simulation.minimizeEnergy()
tetra_nucl.add_reporters(report_interval=10000, output_dcd='output.dcd')
tetra_nucl.simulation.context.setVelocitiesToTemperature(temperature)
tetra_nucl.simulation.step(100000)


