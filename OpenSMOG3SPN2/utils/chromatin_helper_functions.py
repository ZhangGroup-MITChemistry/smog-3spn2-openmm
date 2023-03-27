import numpy as np
import pandas as pd
import sys
import os

'''
Tools for setting up nucleosome related simulations. 
For example, we need functions to remove dihedrals of disordered tails. 
'''

# histone tails (atom index starts from 1): 1-43, 136-159, 238-257, 353-400, 488-530, 623-646, 725-744, 840-887
# in openmm, we use index starting from 0
_histone_tail_start_atoms = np.array([1, 136, 238, 353, 488, 623, 725, 840]) - 1
_histone_tail_end_atoms = np.array([43, 159, 257, 400, 530, 646, 744, 887]) - 1
_histone_tail_atoms = []
for i in range(8):
    _histone_tail_atoms += list(range(_histone_tail_start_atoms[i], _histone_tail_end_atoms[i] + 1))
_histone_tail_atoms = np.array(_histone_tail_atoms)

_n_CA_atoms_per_histone = 974


def remove_histone_tail_dihedrals(df_dihedrals):
    '''
    Remove histone tail dihedral potentials within histone tails. 
    A dihedral potential is removed if at least one atom involved is within histone tail. 
    
    Parameters
    ----------
    df_dihedrals : pd.DataFrame
        Dihedral potential. 
        This should only include histones. 
        Each histone includes 974 CA atoms, and there should be 974*n CG atoms in all (n is the number of histones). 
    
    Returns
    -------
    new_df_dihedrals : pd.DataFrame
        New dihedral potential.
    
    '''
    new_df_dihedrals = pd.DataFrame(columns=df_dihedrals.columns)
    for i, row in df_dihedrals.iterrows():
        a1 = int(row['a1']) % _n_CA_atoms_per_histone
        a2 = int(row['a2']) % _n_CA_atoms_per_histone
        a3 = int(row['a3']) % _n_CA_atoms_per_histone
        a4 = int(row['a4']) % _n_CA_atoms_per_histone
        if not any(x in _histone_tail_atoms for x in [a1, a2, a3, a4]):
            new_df_dihedrals.loc[len(new_df_dihedrals.index)] = row
    return new_df_dihedrals



