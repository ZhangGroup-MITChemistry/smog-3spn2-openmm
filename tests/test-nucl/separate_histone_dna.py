import numpy as np
import pandas as pd
import shutil
import sys
import os
import mdtraj

with open('pdb-files/1kx5.pdb', 'r') as f:
    nucl_pdb_lines = f.readlines()

dna_pdb_lines = []
histone_pdb_lines = []
for each_line in nucl_pdb_lines:
    if each_line[:4] == 'ATOM':
        chainID = each_line[21]
        if chainID in ['I', 'J']:
            dna_pdb_lines.append(each_line)
        elif chainID in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']:
            histone_pdb_lines.append(each_line)

# for dna pdb, clean residue sequence and serial number
r = 1
s = 1
with open('pdb-files/dna.pdb', 'w') as f:
    for i in range(len(dna_pdb_lines)):
        curr_line = dna_pdb_lines[i]
        curr_chainID = curr_line[21]
        curr_resSeq = int(curr_line[22:26])
        if i >= 1:
            prev_line = dna_pdb_lines[i - 1]
            prev_chainID = prev_line[21]
            prev_resSeq = int(prev_line[22:26])
            if curr_resSeq != prev_resSeq:
                r += 1
            if curr_chainID != prev_chainID:
                f.write('TER\n')
                r = 1
        new_line = curr_line[:6] + str(int(s)).rjust(5) + curr_line[11:22] + str(int(r)).rjust(4) + curr_line[26:]
        s += 1
        f.write(new_line)
    f.write('END\n')

s = 1   
with open('pdb-files/histone.pdb', 'w') as f:
    for i in range(len(histone_pdb_lines)):
        curr_line = histone_pdb_lines[i]
        curr_chainID = curr_line[21]
        if i >= 1:
            prev_line = histone_pdb_lines[i - 1]
            prev_chainID = prev_line[21]
            if curr_chainID != prev_chainID:
                f.write('TER\n')
        new_line = curr_line[:6] + str(int(s)).rjust(5) + curr_line[11:]
        s += 1
        f.write(new_line)
    f.write('END\n')


