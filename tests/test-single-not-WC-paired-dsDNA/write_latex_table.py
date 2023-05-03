import numpy as np
import pandas as pd

columns = ['dna bond', 'dna angle', 'dna stacking', 'dna dihedral', 'dna base pair', 'dna cross stacking', 
           'all vdwl', 'all elec']

df1 = pd.read_csv('B_curved_energy.csv')
df2 = pd.read_csv('open3SPN2_B_curved_energy.csv')

for i in range(len(df1.index)):
    line1 = f'    {i + 1} & OpenSMOG3SPN2 & ' + ' & '.join(['%.2f' % x for x in df1.loc[i, columns].tolist()]) + r' \\'
    line2 = f'    {i + 1} & Open3SPN2 & ' + ' & '.join(['%.2f' % x for x in df2.loc[i, columns].tolist()]) + r' \\'
    print(line1)
    print(line2)
    print('    \hline')


