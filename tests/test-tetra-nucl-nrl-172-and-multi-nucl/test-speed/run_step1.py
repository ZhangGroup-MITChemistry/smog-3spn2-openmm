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

'''
Fix tetranucleosome in box center and relax nucleosomes. 

Histone core and middle 73bp of nucleosomal DNA are always rigidified. 
'''

n_single_nucl = 20
box_a, box_b, box_c = 100, 100, 100
platform_name = 'CUDA'
chromatin = SMOG3SPN2Model()

# load tetranucleosome histones
for i in range(4):
    histone_i_parser = SMOGParser.from_atomistic_pdb(f'../tetra-nucl-nrl-172-pdb-files/histone_{i}.pdb', 
                                                     f'histone_{i}_CA.pdb', default_parse=False)
    histone_i_parser.parse_mol(get_native_pairs=False)
    chromatin.append_mol(histone_i_parser)

# load tetranucleosome DNA
with open('tetra_nucl_nrl_172_dna_seq.txt', 'r') as f:
    seq1 = f.readlines()[0].strip()
seq2 = get_WC_paired_seq(seq1)
target_seq = seq1 + seq2

dna_parser = DNA3SPN2Parser.from_atomistic_pdb('../tetra-nucl-nrl-172-pdb-files/dna.pdb', new_sequence=target_seq)
chromatin.append_mol(dna_parser)

n_tetra_nucl_atoms = len(chromatin.atoms.index)
chromatin.atoms_to_pdb('CG_tetra_nucl.pdb')

# load single nucleosomes, first load histone, then load DNA




