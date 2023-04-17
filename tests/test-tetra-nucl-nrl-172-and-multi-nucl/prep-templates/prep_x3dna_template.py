import numpy as np
import pandas as pd
import sys
import os

sys.path.append('../../../..')

from OpenSMOG3SPN2.forcefields.parsers import DNA3SPN2Parser
from OpenSMOG3SPN2.utils.helper_functions import get_WC_paired_seq

'''
Build x3dna templates. So that the templates can be shared with other people, and others can directly load the templates without using x3dna to rebuild. 

'''

# load tetranucleosome DNA
with open('../../tetra_nucl_nrl_172_dna_seq.txt', 'r') as f:
    seq1 = f.readlines()[0].strip()
seq2 = get_WC_paired_seq(seq1)
target_seq = seq1 + seq2

dna_parser = DNA3SPN2Parser.from_atomistic_pdb('../../tetra-nucl-nrl-172-pdb-files/dna.pdb', 
                                               'cg_tetra_nucl_dna.pdb', new_sequence=target_seq, 
                                               temp_name='tetra_nucl_dna')
dna_parser.template_atoms.to_csv('cg_tetra_nucl_dna_template.csv', index=False)

with open('../../single_nucl_dna_seq.txt', 'r') as f:
    seq1 = f.readlines()[0].strip()
seq2 = get_WC_paired_seq(seq1)
target_seq = seq1 + seq2
dna_parser = DNA3SPN2Parser.from_atomistic_pdb('../../single-nucl-pdb-files/dna.pdb', 'cg_single_nucl_dna.pdb', 
                                               new_sequence=target_seq, temp_name='single_nucl_dna')
dna_parser.template_atoms.to_csv('cg_single_nucl_dna_template.csv', index=False)
