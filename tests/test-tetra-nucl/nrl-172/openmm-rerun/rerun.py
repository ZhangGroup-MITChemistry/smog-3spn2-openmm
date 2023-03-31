import numpy as np
import pandas as pd
import sys
import os
import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as unit
import mdtraj

sys.path.append('../../../..')

from OpenSMOG3SPN2.forcefields.parsers import SMOGParser, DNA3SPN2Parser
from OpenSMOG3SPN2.forcefields import SMOG3SPN2Model
from OpenSMOG3SPN2.utils import helper_functions, chromatin_helper_functions

n_nucl = 4
nrl = 172
platform_name = sys.argv[1]
tetra_nucl = SMOG3SPN2Model()
for i in range(n_nucl):
    histone_i_parser = SMOGParser.from_atomistic_pdb(f'../pdb-files/histone_{i}.pdb', f'histone_{i}_CA.pdb', 
                                                     default_parse=False)
    histone_i_parser.parse_mol(get_native_pairs=False)
    old_protein_dihedrals = histone_i_parser.protein_dihedrals
    histone_i_parser.protein_dihedrals = chromatin_helper_functions.remove_histone_tail_dihedrals(old_protein_dihedrals)
    tetra_nucl.append_mol(histone_i_parser)

# load target DNA sequence
# the sequence saved in dna_seq.txt only involves the first ssDNA
# we have to get the full WC-paired sequence
with open('dna_seq.txt', 'r') as f:
    seq1 = f.readlines()[0].strip()
seq2 = helper_functions.get_WC_paired_seq(seq1)
target_seq = seq1 + seq2

dna_parser = DNA3SPN2Parser.from_atomistic_pdb('../pdb-files/dna.pdb', new_sequence=target_seq)
tetra_nucl.append_mol(dna_parser)
tetra_nucl.atoms_to_pdb('cg_chromatin.pdb')

top = app.PDBFile('cg_chromatin.pdb').getTopology()
init_coord = app.PDBFile('cg_chromatin.pdb').getPositions()
tetra_nucl.create_system(top)
tetra_nucl.add_protein_bonds(force_group=1)
tetra_nucl.add_protein_angles(force_group=2)
tetra_nucl.add_protein_dihedrals(force_group=3)
tetra_nucl.add_native_pairs(force_group=4)
tetra_nucl.add_dna_bonds(force_group=5)
tetra_nucl.add_dna_angles(force_group=6)
tetra_nucl.add_dna_stackings(force_group=7)
tetra_nucl.add_dna_dihedrals(force_group=8)
tetra_nucl.add_dna_base_pairs(force_group=9)
tetra_nucl.add_dna_cross_stackings(force_group=10)
tetra_nucl.parse_all_exclusions()
tetra_nucl.exclusions.to_csv('exclusions.csv', index=False)
tetra_nucl.add_all_vdwl(force_group=11)
tetra_nucl.add_all_elec(force_group=12)
tetra_nucl.save_system('system.xml')

temperature = 300*unit.kelvin
#friction_coeff = 0.01/unit.picosecond
timestep = 10*unit.femtosecond
#integrator = mm.LangevinMiddleIntegrator(temperature, friction_coeff, timestep)
integrator = mm.NoseHooverIntegrator(temperature, 1/unit.picosecond, timestep)
tetra_nucl.set_simulation(integrator, platform_name=platform_name, init_coord=init_coord)
simulation = tetra_nucl.simulation

dcd = '../lammps-rerun/traj.dcd'
traj = mdtraj.load_dcd(dcd, top='cg_chromatin.pdb')
n_frames = traj.n_frames
columns = ['protein bond', 'protein angle', 'protein dihedral', 'native pair', 'dna bond', 
           'dna angle', 'dna stacking', 'dna dihedrals', 'dna base pairs', 'dna cross stackings', 
           'all vdwl', 'all electrostatic']
df_energies_kj = pd.DataFrame(columns=columns)
df_energies_kcal = pd.DataFrame(columns=columns)
for i in range(1, n_frames):
    # since lammps gives strange results for the 0th snapshot, we do not compute energy for that one
    simulation.context.setPositions(traj.xyz[i])
    row_kj, row_kcal = [], []
    for j in range(1, 13):
        state = simulation.context.getState(getEnergy=True, groups={j})
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        row_kj.append(energy)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
        row_kcal.append(energy)
        #print(f'Group {j} energy is {energy} kcal/mol')
    df_energies_kj.loc[len(df_energies_kj.index)] = row_kj
    df_energies_kcal.loc[len(df_energies_kcal.index)] = row_kcal

df_energies_kj.round(6).to_csv(f'openmm_energy_kj_{platform_name}.csv', index=False)
df_energies_kcal.round(6).to_csv(f'openmm_energy_kcal_{platform_name}.csv', index=False)


