import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
from scipy.spatial.transform import Rotation as R
import OpenSMOG3SPN2.utils.helper_functions as helper_functions

__author__ = 'Andrew Latham'

__modified_by__ = 'Shuming Liu'

'''
A python script for inserting molecules, similar to gmx insert-molecules. 
'''

def insert_molecules(new_pdb, output_pdb, n_mol, radius=0.5, existing_pdb=None, max_n_attempts=10000, 
                     box=[100, 100, 100], reset_serial=True):
    '''
    Insert multiple copies of given molecule into an existing pdb or a new empty box. 
    Note all the length parameter units are nm, though coordinates in pdb are in unit angstroms. 
    Currently this function only supports using orthogonal box. 
    
    Parameters
    ----------
    new_pdb : str
        Path of the new pdb file to be inserted. 
    
    output_pdb : str
        Path of the output pdb file. 
    
    n_mol : int
        Number of copies of the new molecule to be inserted. 
    
    radius : float or int
        Radius for each atom in unit nm. Insertion ensures the distance between any two atoms are larger than 2*radius. 
    
    existing_pdb : None or str
        Path of an existing pdb file. If None, the molecules will be inserted into a new empty box, otherwise, inserted into a copy of existing_pdb. 
    
    max_n_attempts : int
        Maximal number of attempts to insert. Ensure this value no smaller than n_mol. 
    
    box : list or numpy array, shape = (3,)
        New box size in unit nm. This is effective if reset_box is True. 
    
    reset_serial : bool
        Whether to reset serial to 0, 1, ..., N - 1. 
    
    '''
    if existing_pdb is None:
        atoms = pd.DataFrame()
    else:
        atoms = helper_functions.parse_pdb(existing_pdb)
    new_atoms = helper_functions.parse_pdb(new_pdb)
    new_coord = new_atoms[['x', 'y', 'z']].to_numpy()
    new_coord -= np.mean(new_coord, axis=0)
    assert max_n_attempts >= n_mol
    count_n_mol = 0
    count_n_attempts = 0
    cutoff = float(2*10*radius)
    box_a, box_b, box_c = 10*box[0], 10*box[1], 10*box[2] # convert nm to angstrom
    dim = np.array([box_a, box_b, box_c, 90.0, 90.0, 90.0])
    while (count_n_mol < n_mol) and (count_n_attempts < max_n_attempts):
        # get a random rotation
        rotate = R.random()
        new_coord_i = rotate.apply(new_coord)
        # get a random translation
        translate = np.random.uniform(0, 1, 3)*np.array([box_a, box_b, box_c])
        new_coord_i += translate
        if len(atoms.index) == 0:
            new_atoms_i = new_atoms.copy()
            new_atoms_i[['x', 'y', 'z']] = new_coord_i
            atoms = pd.concat([atoms, new_atoms_i], ignore_index=True)
            count_n_mol += 1
        else:
            coord = atoms[['x', 'y', 'z']].to_numpy()
            d = distance_array(coord, new_coord_i, dim)
            if np.amin(d) >= cutoff:
                new_atoms_i = new_atoms.copy()
                new_atoms_i[['x', 'y', 'z']] = new_coord_i
                atoms = pd.concat([atoms, new_atoms_i], ignore_index=True)
                count_n_mol += 1
        count_n_attempts += 1
    
    # determine if the number of molecules were successfully added:
    if count_n_mol == n_mol:
        print(f'Successfully inserted {n_mol} molecules.')
    else:
        print(f'Could not successfully insert {n_mol} molecules in {count_n_attempts} attempts.')
        print(f'Only added {count_n_mol} molecules. Try increasing the box size or number of attempts to add more molecules.')
    if reset_serial:
        atoms['serial'] = list(range(len(atoms.index)))
    # write the final pdb
    helper_functions.write_pdb(atoms, output_pdb)


def insert_molecules_beta(new_pdb, output_pdb, n_mol, radius=0.5, existing_pdb=None, max_n_attempts=10000, 
                          box=[100, 100, 100], reset_serial=True):
    '''
    This is a beta version.
    This should be faster then `insert_molecules` for large systems. 
    
    Insert multiple copies of given molecule into an existing pdb or a new empty box. 
    Note all the length parameter units are nm, though coordinates in pdb are in unit angstroms. 
    Currently this function only supports using orthogonal box. 
    
    Parameters
    ----------
    new_pdb : str
        Path of the new pdb file to be inserted. 
    
    output_pdb : str
        Path of the output pdb file. 
    
    n_mol : int
        Number of copies of the new molecule to be inserted. 
    
    radius : float or int
        Radius for each atom in unit nm. Insertion ensures the distance between any two atoms are larger than 2*radius. 
    
    existing_pdb : None or str
        Path of an existing pdb file. If None, the molecules will be inserted into a new empty box, otherwise, inserted into a copy of existing_pdb. 
    
    max_n_attempts : int
        Maximal number of attempts to insert. Ensure this value no smaller than n_mol. 
    
    box : list or numpy array, shape = (3,)
        New box size in unit nm. This is effective if reset_box is True. 
    
    reset_serial : bool
        Whether to reset serial to 0, 1, ..., N - 1. 
    
    '''
    if existing_pdb is None:
        atoms = pd.DataFrame()
    else:
        atoms = helper_functions.parse_pdb(existing_pdb)
    new_atoms = helper_functions.parse_pdb(new_pdb)
    new_coord = new_atoms[['x', 'y', 'z']].to_numpy()
    new_coord -= np.mean(new_coord, axis=0)
    assert max_n_attempts >= n_mol
    count_n_mol = 0
    count_n_attempts = 0
    cutoff = float(2*10*radius)
    box_a, box_b, box_c = 10*box[0], 10*box[1], 10*box[2] # convert nm to angstrom
    dim = np.array([box_a, box_b, box_c, 90.0, 90.0, 90.0]).astype(np.float32)
    if len(atoms.index) > 0:
        coord = atoms[['x', 'y', 'z']].to_numpy().astype(np.float32)
        grid_search = mda.lib.nsgrid.FastNS(cutoff, coord, dim, pbc=True)
    while (count_n_mol < n_mol) and (count_n_attempts < max_n_attempts):
        # get a random rotation
        rotate = R.random()
        new_coord_i = rotate.apply(new_coord)
        # get a random translation
        translate = np.random.uniform(0, 1, 3)*np.array([box_a, box_b, box_c])
        new_coord_i += translate
        # test with mdanalysis 2.4.3
        # if new_coord_i dtype is np.float64, there will be error
        # if new_coord_i dtype is np.float32, the code can run smoothly
        new_coord_i = new_coord_i.astype(np.float32)
        if len(atoms.index) == 0:
            new_atoms_i = new_atoms.copy()
            new_atoms_i[['x', 'y', 'z']] = new_coord_i
            atoms = pd.concat([atoms, new_atoms_i], ignore_index=True)
            count_n_mol += 1
            coord = atoms[['x', 'y', 'z']].to_numpy().astype(np.float32)
            grid_search = mda.lib.nsgrid.FastNS(cutoff, coord, dim, pbc=True)
        else:
            results = grid_search.search(new_coord_i)
            if len(results.get_pair_distances()) == 0:
                new_atoms_i = new_atoms.copy()
                new_atoms_i[['x', 'y', 'z']] = new_coord_i
                atoms = pd.concat([atoms, new_atoms_i], ignore_index=True)
                count_n_mol += 1
                coord = atoms[['x', 'y', 'z']].to_numpy().astype(np.float32)
                grid_search = mda.lib.nsgrid.FastNS(cutoff, coord, dim, pbc=True)
        count_n_attempts += 1
    
    # determine if the number of molecules were successfully added:
    if count_n_mol == n_mol:
        print(f'Successfully inserted {n_mol} molecules.')
    else:
        print(f'Could not successfully insert {n_mol} molecules in {count_n_attempts} attempts.')
        print(f'Only added {count_n_mol} molecules. Try increasing the box size or number of attempts to add more molecules.')
    if reset_serial:
        atoms['serial'] = list(range(len(atoms.index)))
    # write the final pdb
    helper_functions.write_pdb(atoms, output_pdb)
    
