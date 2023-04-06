# %%
import numpy as np
import pandas as pd
import sys
import os
import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as unit
import mdtraj

sys.path.append('../..')

from OpenSMOG3SPN2.forcefields.parsers import SMOGParser, DNA3SPN2Parser
from OpenSMOG3SPN2.forcefields import SMOG3SPN2Model
from OpenSMOG3SPN2.utils.helper_functions import get_WC_paired_seq
from OpenSMOG3SPN2.utils.chromatin_helper_functions import remove_histone_tail_dihedrals, get_chromatin_rigid_bodies
from OpenSMOG3SPN2.forcefields.rigid import createRigidBodies

# %% [markdown]
# This is the tutorial for setting tetranucleosome simulation. This tetranucleosome has NRL value equal to 172. 

# %%
# set parameters
n_nucl = 4
nrl = 172

# load histones
# pdb-files/histone_{i}.pdb is the pdb of the i-th histone
# as histone core is rigid, we do not need to search native pairs
tetra_nucl = SMOG3SPN2Model()
for i in range(n_nucl):
    histone_i_parser = SMOGParser.from_atomistic_pdb(f'pdb-files/histone_{i}.pdb', f'histone_{i}_CA.pdb', 
                                                     default_parse=False)
    histone_i_parser.parse_mol(get_native_pairs=False)
    # remove dihedrals involving histone tails atoms
    histone_i_parser.protein_dihedrals = remove_histone_tail_dihedrals(histone_i_parser.protein_dihedrals)
    tetra_nucl.append_mol(histone_i_parser)

# load dna with target sequence
# note for DNA, full sequence is required, as 3SPN2 does not require dsDNA seuqence to be W-C paired
with open('dna_seq.txt', 'r') as f:
    seq1 = f.readlines()[0].strip()
seq2 = get_WC_paired_seq(seq1)
target_seq = seq1 + seq2

dna_parser = DNA3SPN2Parser.from_atomistic_pdb('pdb-files/dna.pdb', new_sequence=target_seq)
tetra_nucl.append_mol(dna_parser)
tetra_nucl.atoms_to_pdb('cg_chromatin.pdb')


# %% [markdown]
# Prepare openmm system, add forces, and run simulation. 

# %%
# prepare system
# set rigid body with the coordinates of the final snapshot of traj.dcd
top = app.PDBFile('cg_chromatin.pdb').getTopology()
traj = mdtraj.load_dcd('traj.dcd', 'cg_chromatin.pdb')
init_coord = traj.xyz[-1]*unit.nanometer
rigid_coord = init_coord
tetra_nucl.create_system(top, box_a=200, box_b=200, box_c=200)

# parse exclusions and set rigid bodies
# though we parse exclusions before setting rigid bodies, we do not change exclusion list when we set rigid bodies, thus avoiding bonded atoms within rigid bodies have extremely strong nonbonded interactions
tetra_nucl.parse_all_exclusions()
rigid_bodies = get_chromatin_rigid_bodies(n_nucl, nrl, n_rigid_bp_per_nucl=73)
tetra_nucl.set_rigid_bodies(rigid_coord, rigid_bodies, keep_unchanged=['protein_exclusions', 'dna_exclusions', 
                                                                       'exclusions'])


# %%
# add forces
tetra_nucl.add_protein_bonds(force_group=1)
tetra_nucl.add_protein_angles(force_group=2)
#tetra_nucl.add_protein_dihedrals(force_group=3) # no need to add dihedrals as they are all within the rigid bodies
#tetra_nucl.add_native_pairs(force_group=4) # no need to add native pairs as they are all within the rigid bodies
tetra_nucl.add_dna_bonds(force_group=5)
tetra_nucl.add_dna_angles(force_group=6)
tetra_nucl.add_dna_stackings(force_group=7)
tetra_nucl.add_dna_dihedrals(force_group=8)
tetra_nucl.add_dna_base_pairs(force_group=9)
tetra_nucl.add_dna_cross_stackings(force_group=10)
tetra_nucl.add_all_vdwl(force_group=11)
tetra_nucl.add_all_elec(force_group=12)
tetra_nucl.save_system('rigid_system.xml')

# %%
# set and run simulation
temperature = 300*unit.kelvin
friction_coeff = 0.01/unit.picosecond
timestep = 5*unit.femtosecond
integrator = mm.LangevinMiddleIntegrator(temperature, friction_coeff, timestep)
tetra_nucl.set_simulation(integrator, platform_name='CUDA', init_coord=init_coord)

print('Simulation starting configruation energies of each group:')
for i in range(1, 13):
    state = tetra_nucl.simulation.context.getState(getEnergy=True, groups={i})
    energy = state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
    print(f'Group {i} energy is {energy} kcal/mol.')

tetra_nucl.simulation.minimizeEnergy()
tetra_nucl.add_reporters(report_interval=100, output_dcd='output.dcd')
tetra_nucl.simulation.context.setVelocitiesToTemperature(temperature)
tetra_nucl.simulation.step(1000)


