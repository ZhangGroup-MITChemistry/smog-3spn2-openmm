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
from OpenSMOG3SPN2.utils.helper_functions import get_WC_paired_seq, write_pdb
from OpenSMOG3SPN2.utils.chromatin_helper_functions import get_chromatin_rigid_bodies
from OpenSMOG3SPN2.utils.insert import insert_molecules

'''
Build system xml file without any rigid bodies.

We remove bonds, angles, protein 1-3 and 1-4 exclusions within rigid bodies, so that the system is smaller. 
'''

if not os.path.exists('intermediate-files'):
    os.makedirs('intermediate-files')

n_single_nucl = 26
box_a, box_b, box_c = 55.0, 55.0, 55.0
chromatin = SMOG3SPN2Model()

# load tetranucleosome histones
for i in range(4):
    histone_i_parser = SMOGParser.from_atomistic_pdb(f'../tetra-nucl-nrl-172-pdb-files/histone_{i}.pdb', 
                                                     f'intermediate-files/cg_tetra_nucl_histone_{i}.pdb', 
                                                     default_parse=False)
    histone_i_parser.parse_mol(get_native_pairs=False)
    chromatin.append_mol(histone_i_parser)

# load tetranucleosome DNA
with open('tetra_nucl_nrl_172_dna_seq.txt', 'r') as f:
    seq1 = f.readlines()[0].strip()
seq2 = get_WC_paired_seq(seq1)
target_seq = seq1 + seq2

dna_parser = DNA3SPN2Parser.from_atomistic_pdb('../tetra-nucl-nrl-172-pdb-files/dna.pdb', 
                                               'intermediate-files/cg_tetra_nucl_dna.pdb', new_sequence=target_seq, 
                                               temp_name='temp1')
chromatin.append_mol(dna_parser)

chromatin.atoms_to_pdb('intermediate-files/cg_tetra_nucl_original_coord.pdb')
n_tetra_nucl_atoms = len(chromatin.atoms.index)
tetra_nucl_atoms = chromatin.atoms.copy()
# put tetranucleosome geometric center at the center of the box
traj = mdtraj.load_dcd('../tetra_nucl_nrl_172_traj.dcd', top='intermediate-files/cg_tetra_nucl_original_coord.pdb')
coord = traj.xyz[-1]
coord -= np.mean(coord, axis=0)
coord += np.array([box_a, box_b, box_c])/2
tetra_nucl_atoms[['x', 'y', 'z']] = 10*coord # convert nm to A
tetra_nucl_atoms.loc[:, 'charge'] = ''
write_pdb(tetra_nucl_atoms, 'cg_tetra_nucl.pdb')

# load single nucleosomes, first load histone, then load DNA
histone_parser = SMOGParser.from_atomistic_pdb('../single-nucl-pdb-files/histone.pdb', 
                                               'intermediate-files/cg_single_nucl_histone.pdb', default_parse=False)
histone_parser.parse_mol(get_native_pairs=False)
with open('single_nucl_dna_seq.txt', 'r') as f:
    seq1 = f.readlines()[0].strip()
seq2 = get_WC_paired_seq(seq1)
target_seq = seq1 + seq2
dna_parser = DNA3SPN2Parser.from_atomistic_pdb('../single-nucl-pdb-files/dna.pdb', 
                                               'intermediate-files/cg_single_nucl_dna.pdb', new_sequence=target_seq, 
                                               temp_name='temp2')
nucl = SMOG3SPN2Model()
nucl.append_mol(histone_parser)
nucl.append_mol(dna_parser)
nucl.atoms_to_pdb('cg_single_nucl.pdb')

# prepare pdb file for the whole system
if not os.path.exists('start.pdb'):
    insert_molecules('cg_single_nucl.pdb', 'start.pdb', n_mol=n_single_nucl, existing_pdb='cg_tetra_nucl.pdb', 
                     box=[box_a, box_b, box_c])
top = app.PDBFile('start.pdb').getTopology()

for i in range(n_single_nucl):
    chromatin.append_mol(histone_parser)
    chromatin.append_mol(dna_parser)

# set rigid bodies
# rigid bodies are histone cores with middle 73 bp of core DNA
# we remove bonds, angles, and 1-3 exclusions within the rigid bodies
rigid_bodies = []
n_tetra_nucl_atoms = app.PDBFile('cg_tetra_nucl.pdb').getTopology().getNumAtoms()
tetra_nucl_rigid_bodies = get_chromatin_rigid_bodies(n_nucl=4, nrl=172)
rigid_bodies += tetra_nucl_rigid_bodies
n_atoms_per_single_nucl = app.PDBFile('cg_single_nucl.pdb').getTopology().getNumAtoms()
single_nucl_rigid_body = np.array(get_chromatin_rigid_bodies(n_nucl=1, nrl=147)[0])
for i in range(n_single_nucl):
    rigid_bodies.append((single_nucl_rigid_body + n_tetra_nucl_atoms + i*n_atoms_per_single_nucl).tolist())

n_atoms = app.PDBFile('start.pdb').getTopology().getNumAtoms()
rigid_body_dict = {}
for i in range(n_atoms):
    rigid_body_dict[i] = None
for i in range(len(rigid_bodies)):
    for j in rigid_bodies[i]:
        rigid_body_dict[j] = i

new_protein_bonds = pd.DataFrame(columns=chromatin.protein_bonds.columns)
for i, row in chromatin.protein_bonds.iterrows():
    a1, a2 = int(row['a1']), int(row['a2'])
    r1, r2 = rigid_body_dict[a1], rigid_body_dict[a2]
    if (r1 is not None) and (r1 == r2):
        pass
    else:
        new_protein_bonds.loc[len(new_protein_bonds.index)] = row
chromatin.protein_bonds = new_protein_bonds

new_protein_angles = pd.DataFrame(columns=chromatin.protein_angles.columns)
for i, row in chromatin.protein_angles.iterrows():
    a1, a2, a3 = int(row['a1']), int(row['a2']), int(row['a3'])
    r1, r2, r3 = rigid_body_dict[a1], rigid_body_dict[a2], rigid_body_dict[a3]
    if (r1 is not None) and (r1 == r2) and (r1 == r3):
        pass
    else:
        new_protein_angles.loc[len(new_protein_angles.index)] = row
chromatin.protein_angles = new_protein_angles

new_protein_exclusions = pd.DataFrame(columns=chromatin.protein_exclusions.columns)
for i, row in chromatin.protein_exclusions.iterrows():
    a1, a2 = int(row['a1']), int(row['a2'])
    r1, r2 = rigid_body_dict[a1], rigid_body_dict[a2]
    if (abs(a1 - a2) >= 2) and (r1 is not None) and (r1 == r2):
        pass
    else:
        new_protein_exclusions.loc[len(new_protein_exclusions.index)] = row
chromatin.protein_exclusions = new_protein_exclusions

chromatin.create_system(top, box_a=box_a, box_b=box_b, box_c=box_c)
chromatin.add_protein_bonds(force_group=1)
chromatin.add_protein_angles(force_group=2)
chromatin.add_dna_bonds(force_group=5)
chromatin.add_dna_angles(force_group=6)
chromatin.add_dna_stackings(force_group=7)
chromatin.add_dna_dihedrals(force_group=8)
chromatin.add_dna_base_pairs(force_group=9)
chromatin.add_dna_cross_stackings(force_group=10)
chromatin.parse_all_exclusions()
chromatin.add_all_vdwl(force_group=11)
chromatin.add_all_elec(force_group=12)
chromatin.save_system('nonrigid_light_system.xml')


