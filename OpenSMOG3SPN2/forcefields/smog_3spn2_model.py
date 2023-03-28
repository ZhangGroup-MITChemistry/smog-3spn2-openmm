import numpy as np
import pandas as pd
import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as unit
import itertools
from OpenSMOG3SPN2.forcefields.cg_model import CGModel
from OpenSMOG3SPN2.forcefields import functional_terms
from OpenSMOG3SPN2.forcefields.parameters.mixin_3spn2_config_parser import Mixin3SPN2ConfigParser
import sys
import os

__location__ = os.path.dirname(os.path.abspath(__file__))

_amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS',
                'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
                'LEU', 'LYS', 'MET', 'PHE', 'PRO',
                'SER', 'THR', 'TRP', 'TYR', 'VAL']

_nucleotides = ['DA', 'DT', 'DC', 'DG']

_WC_pair_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

class SMOG3SPN2Model(CGModel, Mixin3SPN2ConfigParser):
    '''
    A class for SMOG+3SPN2 model. 
    To ensure this model works properly, please ensure two neighboring ssDNA chains do not share same chainID. 
    '''
    def __init__(self, dna_type='B_curved', default_parse_config=True):
        '''
        Initialize. 
        
        Parameters
        ----------
        default_parse_config : bool
            Whether to parse the default 3SPN2 configuration file. 
        
        dna_type : str
            DNA type.
        
        '''
        self.atoms = None
        # note we set base pair and cross stacking donors and acceptors as bonded attributes
        self.bonded_attr_names = ['protein_bonds', 'protein_angles', 'protein_dihedrals', 'native_pairs', 
                                  'dna_bonds', 'dna_angles', 'dna_stackings', 'dna_dihedrals', 'base_pair_donor_A', 
                                  'base_pair_donor_G', 'base_pair_acceptor_C', 'base_pair_acceptor_T', 'exclusions']
        self.dna_type = dna_type
        if default_parse_config:
            # load parameters
            self.parse_config_file()
    
    # rewrite append_mol method as we need to check DNA type consistency when appending new DNA
    def append_mol(self, new_mol, verbose=False):
        '''
        The method can append new molecules by concatenating atoms and bonded interaction information saved in dataframes. 
        Ensure two neighboring chains do not share chainID. 
        
        Parameters
        ----------
        new_mol : a consistent parser object or CG model object
            The object of a new molecule including atom and bonded interaction information. 
        
        verbose : bool
            Whether to report the appended attributes. 
        
        '''
        if hasattr(new_mol, 'dna_type'):
            # set DNA type
            if hasattr(self, 'dna_type'):
                if getattr(self, 'dna_type') is None:
                    setattr(self, 'dna_type', getattr(new_mol, 'dna_type'))
                else:
                    # check DNA type consistency
                    assert getattr(self, 'dna_type') == getattr(new_mol, 'dna_type')
            else:
                setattr(self, 'dna_type', getattr(new_mol, 'dna_type'))
        
        new_atoms = new_mol.atoms.copy()
        if hasattr(self, 'atoms'):
            if self.atoms is None:
                add_index = 0
                self.atoms = new_atoms
            else:
                add_index = len(self.atoms.index)
                self.atoms = pd.concat([self.atoms, new_atoms], ignore_index=True)
        else:
            add_index = 0
            self.atoms = new_atoms
        for each_attr_name in self.bonded_attr_names:
            if verbose:
                print(f'Append attribute: {each_attr_name}. ')
            if hasattr(new_mol, each_attr_name):
                if getattr(new_mol, each_attr_name) is not None:
                    new_attr = getattr(new_mol, each_attr_name).copy()
                    for each_col in ['a1', 'a2', 'a3', 'a4']:
                        if each_col in new_attr.columns:
                            new_attr[each_col] += add_index
                    if hasattr(self, each_attr_name):
                        if getattr(self, each_attr_name) is None:
                            setattr(self, each_attr_name, new_attr)
                        else:
                            combined_attr = pd.concat([getattr(self, each_attr_name).copy(), new_attr], 
                                                      ignore_index=True)
                            setattr(self, each_attr_name, combined_attr)
                    else:
                        setattr(self, each_attr_name, new_attr)

    def add_protein_bonds(self, force_group=1):
        '''
        Add protein bonds.
        
        Parameters
        ----------
        force_group : int
            Force group. 
        
        '''
        if hasattr(self, 'protein_bonds'):
            print('Add protein bonds.')
            force = functional_terms.harmonic_bond_term(self.protein_bonds, self.use_pbc, force_group)
            self.system.addForce(force)
    
    def add_protein_angles(self, force_group=2):
        '''
        Add protein angles.
        
        Parameters
        ----------
        force_group : int
            Force group. 
        
        '''
        if hasattr(self, 'protein_angles'):
            print('Add protein angles.')
            force = functional_terms.harmonic_angle_term(self.protein_angles, self.use_pbc, force_group)
            self.system.addForce(force)

    def add_protein_dihedrals(self, force_group=3):
        '''
        Add protein dihedrals. 
        
        Parameters
        ----------
        force_group : int
            Force group. 
        
        '''
        if hasattr(self, 'protein_dihedrals'):
            print('Add protein dihedrals.')
            force = functional_terms.periodic_dihedral_term(self.protein_dihedrals, self.use_pbc, force_group)
            self.system.addForce(force)

    def add_native_pairs(self, force_group=4):
        '''
        Add native pairs. 
        
        Parameters
        ----------
        force_group : int
            Force group.
        
        '''
        if hasattr(self, 'native_pairs'):
            print('Add native pairs.')
            force = functional_terms.native_pair_gaussian_term(self.native_pairs, self.use_pbc, force_group)
            self.system.addForce(force)
    
    def add_dna_bonds(self, force_group=5):
        '''
        Add DNA bonds.
        
        Parameters
        ----------
        force_group : int
            Force group.
        
        '''
        if hasattr(self, 'dna_bonds'):
            print('Add DNA bonds.')
            force = functional_terms.class2_bond_term(self.dna_bonds, self.use_pbc, force_group)
            self.system.addForce(force)
        
    def add_dna_angles(self, force_group=6):
        '''
        Add DNA angles. 
        
        Parameters
        ----------
        force_group : int
            Force group.
        
        '''
        if hasattr(self, 'dna_angles'):
            print('Add DNA angles.')
            force = functional_terms.harmonic_angle_term(self.dna_angles, self.use_pbc, force_group)
            self.system.addForce(force)
        
    def add_dna_stackings(self, force_group=7):
        '''
        Add DNA stackings. 
        
        Parameters
        ----------
        force_group : int
            Force group.
        
        '''
        if hasattr(self, 'dna_stackings'):
            print('Add DNA stackings.')
            force = functional_terms.dna_3spn2_stacking_term(self.dna_stackings, self.use_pbc, force_group)
            self.system.addForce(force)
    
    def add_dna_dihedrals(self, force_group=8):
        '''
        Add DNA dihedrals.
        
        Parameters
        ----------
        force_group : int
            Force group.
        
        '''
        if hasattr(self, 'dna_dihedrals'):
            print('Add DNA dihedrals.')
            force = functional_terms.dna_3spn2_dihedral_term(self.dna_dihedrals, self.use_pbc, force_group)
            self.system.addForce(force)
    
    def add_dna_base_pairs(self, cutoff=1.8, force_group=9):
        '''
        Add DNA base pair potentials. 
        
        Parameters
        ----------
        cutoff : float or int
            Cutoff distance. 
        
        force_group : int
            Force group.
        
        '''
        print('Add DNA base pairs.')
        pair_definition = self.pair_definition[self.pair_definition['DNA'] == self.dna_type]
        atoms = self.atoms.copy()
        # reset chainID to unique numbers so forces can be set properly
        # please ensure in input self.atoms two neighboring chains do not share the same chainID
        # this is because the method recognizes individual chains by comparing chainIDs of neighboring atoms
        atoms.index = list(range(len(atoms.index)))
        new_chainIDs = []
        c = 0
        for i, row in atoms.iterrows():
            if i >= 1:
                if row['chainID'] != atoms.loc[i - 1, 'chainID']:
                    c += 1
            new_chainIDs.append(c)
        atoms['chainID'] = new_chainIDs
        # change index so we can locate each atom easily
        atoms['index'] = list(range(len(atoms.index)))
        atoms.index = zip(atoms['chainID'], atoms['resSeq'], atoms['name'])
        dna_atoms = atoms[atoms['resname'].isin(_nucleotides)]
        # find W-C pairs
        # note as we have two types of W-C pairs (A-T and C-G), we need two separate forces
        for i, row in pair_definition.iterrows():
            parameters = [row['torsion'], row['sigma'], row['t1'], row['t2'], row['rang'], 
                          row['epsilon'], row['alpha']]
            base1, base2 = row['Base1'], row['Base2']
            donors1 = dna_atoms[dna_atoms['name'] == base1].copy()
            acceptors1 = dna_atoms[dna_atoms['name'] == base2].copy()
            donors2 = dna_atoms.loc[[(x[0], x[1], 'S') for x in donors1.index]]
            acceptors2 = dna_atoms.loc[[(x[0], x[1], 'S') for x in acceptors1.index]]
            force_i = functional_terms.dna_3spn2_base_pair_term(self.use_pbc, cutoff, force_group)
            # add donors and acceptors
            for a1, a2 in zip(donors1['index'].tolist(), donors2['index'].tolist()):
                force_i.addDonor(a1, a2, -1, parameters)
            for a1, a2 in zip(acceptors1['index'].tolist(), acceptors2['index'].tolist()):
                force_i.addAcceptor(a1, a2, -1)
            # set exclusions
            # for hbond forces, exclusions are added based on donor id and acceptor id
            donors1['donor_id'] = list(range(len(donors1.index)))
            acceptors1['acceptor_id'] = list(range(len(acceptors1.index)))
            dna_unique_chainIDs = dna_atoms['chainID'].drop_duplicates(keep='first').tolist()
            for c in dna_unique_chainIDs:
                donors1_c = donors1[donors1['chainID'] == c]
                acceptors1_c = acceptors1[acceptors1['chainID'] == c]
                for j, atom1 in donors1_c.iterrows():
                    for k, atom2 in acceptors1_c.iterrows():
                        if abs(atom1['resSeq'] - atom2['resSeq']) <= 2:
                            force_i.addExclusion(atom1['donor_id'], atom2['acceptor_id'])
            self.system.addForce(force_i)
    
    def add_dna_cross_stackings(self, cutoff=1.8, OpenCLPatch=True, force_group=10):
        '''
        Add DNA cross stacking potentials. 
        '''
        print('Add DNA cross stackings.')
        cross_definition = self.cross_definition[self.cross_definition['DNA'] == self.dna_type].copy()
        cross_definition.index = [x for x in zip(cross_definition['Base_d1'], cross_definition['Base_a1'], 
                                                 cross_definition['Base_a3'])]
        atoms = self.atoms.copy()
        # reset chainID to unique numbers so forces can be set properly
        # please ensure in input self.atoms two neighboring chains do not share the same chainID
        # this is because the method recognizes individual chains by comparing chainIDs of neighboring atoms
        atoms.index = list(range(len(atoms.index)))
        new_chainIDs = []
        c = 0
        for i, row in atoms.iterrows():
            if i >= 1:
                if row['chainID'] != atoms.loc[i - 1, 'chainID']:
                    c += 1
            new_chainIDs.append(c)
        atoms['chainID'] = new_chainIDs
        # change index so we can locate each atom easily
        atoms['index'] = list(range(len(atoms.index)))
        dna_atoms = atoms[atoms['resname'].isin(_nucleotides)].copy()
        dna_atoms['group'] = dna_atoms['name']
        dna_atoms.loc[dna_atoms['name'].isin(['A', 'T', 'C', 'G']), 'group'] = 'B'
        dna_atoms.index = zip(dna_atoms['chainID'], dna_atoms['resSeq'], dna_atoms['group'])
        bases = dna_atoms[dna_atoms['group'] == 'B']
        sugars = dna_atoms.loc[[(x[0], x[1], 'S') for x in bases.index]]
        next_bases = bases.reindex([(x[0], x[1] + 1, 'B') for x in bases.index])
        prev_bases = bases.reindex([(x[0], x[1] - 1, 'B') for x in bases.index])
        bases.index = list(range(len(bases.index)))
        sugars.index = bases.index
        next_bases.index = bases.index
        prev_bases.index = bases.index
        sel_next_bases = next_bases[next_bases['name'].isin(['A', 'T', 'C', 'G'])].index.tolist()
        sel_prev_bases = prev_bases[prev_bases['name'].isin(['A', 'T', 'C', 'G'])].index.tolist()
        # define forces
        dict_cross_stackings = {}
        for b in ['A', 'T', 'G', 'C']:
            force1 = functional_terms.dna_3spn2_cross_stacking_term(self.use_pbc, cutoff, force_group)
            force2 = functional_terms.dna_3spn2_cross_stacking_term(self.use_pbc, cutoff, force_group)
            dict_cross_stackings.update({b: (force1, force2)})
        donor_dict = {i: [] for i in ['A', 'T', 'G', 'C']}
        acceptor_dict = {i: [] for i in ['A', 'T', 'G', 'C']}
        # use a1, a2, a3 for atom indices
        for i in sel_next_bases:
            a1, a1n = bases.loc[i, 'index'], bases.loc[i, 'name']
            a2 = sugars.loc[i, 'index']
            a3, a3n = next_bases.loc[i, 'index'], next_bases.loc[i, 'name']
            a1, a2, a3 = int(a1), int(a2), int(a3)
            force1, force2 = dict_cross_stackings[a1n]
            p = cross_definition.loc[[(_WC_pair_dict[a1n], a1n, a3n)]].squeeze()
            p = p[['t03', 'T0CS_2', 'rng_cs2', 'rng_bp', 'eps_cs2', 'alpha_cs2', 'Sigma_2']].tolist()
            force1.addDonor(a1, a2, a3)
            force2.addAcceptor(a1, a2, a3, p)
            donor_dict[a1n].append(a1) # add force1 donor a1
        for i in sel_prev_bases:
            a1, a1n = bases.loc[i, 'index'], bases.loc[i, 'name']
            a2 = sugars.loc[i, 'index']
            a3, a3n = prev_bases.loc[i, 'index'], prev_bases.loc[i, 'name']
            a1, a2, a3 = int(a1), int(a2), int(a3)
            force1, force2 = dict_cross_stackings[_WC_pair_dict[a1n]]
            p = cross_definition.loc[[(_WC_pair_dict[a1n], a1n, a3n)]].squeeze()
            p = p[['t03', 'T0CS_1', 'rng_cs1', 'rng_bp', 'eps_cs1', 'alpha_cs1', 'Sigma_1']].tolist()
            force1.addAcceptor(a1, a2, a3, p)
            force2.addDonor(a1, a2, a3)
            acceptor_dict[_WC_pair_dict[a1n]].append(a1) # add force1 acceptor a1
        # set exclusions
        dna_unique_chainIDs = dna_atoms['chainID'].drop_duplicates(keep='first').tolist()
        if OpenCLPatch:
            max_n = 6
        else:
            max_n = 9
        for b in ['A', 'T', 'C', 'G']:
            force1, force2 = dict_cross_stackings[b]
            donors = dna_atoms[dna_atoms['index'].isin(donor_dict[b])].copy()
            acceptors = dna_atoms[dna_atoms['index'].isin(acceptor_dict[b])].copy()
            donors['donor_id'] = list(range(len(donors.index)))
            acceptors['acceptor_id'] = list(range(len(acceptors.index)))
            for c in dna_unique_chainIDs:
                donors_c = donors[donors['chainID'] == c]
                acceptors_c = acceptors[acceptors['chainID'] == c]
                for j, atom1 in donors_c.iterrows():
                    for k, atom2 in acceptors_c.iterrows():
                        a1, a2 = atom1['index'], atom2['index']
                        if (abs(a1 - a2) <= max_n) or ((not OpenCLPatch) and (a1 > a2)):
                            force1.addExclusion(atom1['donor_id'], atom2['acceptor_id'])
                            force2.addExclusion(atom2['acceptor_id'], atom1['donor_id'])
            self.system.addForce(force1)
            self.system.addForce(force2)

            
                    
            
            
            
        
            
            
            
            
            
            
            
            
        
        
        
        
        

    
