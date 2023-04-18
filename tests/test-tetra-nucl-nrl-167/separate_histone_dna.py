import numpy as np
import shutil
import sys
import os

# separate dsDNA and histones into individual pdb files
n_nucl = 4
dna_pdb_path = 'pdb-files/fiber-167-4_clean.pdb'
histones_pdb_path = 'pdb-files/histone-all.pdb'

shutil.copyfile(dna_pdb_path, 'pdb-files/dna.pdb')

with open(histones_pdb_path, 'r') as f:
    histones_pdb_lines = f.readlines()
n_histones_pdb_lines = len(histones_pdb_lines)

output_histone_pdb_lines = []
for i in range(n_nucl):
    output_histone_pdb_lines.append([])

histone_index = 0
for i in range(n_histones_pdb_lines):
    line_i = histones_pdb_lines[i]
    if line_i[:4] == 'ATOM':
        output_histone_pdb_lines[histone_index].append(line_i)
    if line_i[:3] == 'TER':
        curr_chain_id = histones_pdb_lines[i - 1][21]
        if curr_chain_id == 'H':
            output_histone_pdb_lines[histone_index].append('END\n')
            histone_index += 1
        else:
            output_histone_pdb_lines[histone_index].append(line_i)
    if line_i[:3] == 'END':
        output_histone_pdb_lines[histone_index].append('END\n')
            
for i in range(n_nucl):
    with open(f'pdb-files/histone_{i}.pdb', 'w') as f:
        for each_line in output_histone_pdb_lines[i]:
            f.write(each_line)


