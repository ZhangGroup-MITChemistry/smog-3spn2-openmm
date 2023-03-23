import numpy as np
import pandas as pd
import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as unit
from .cg_model import CGModel
from . import functional_terms
import sys
import os

__location__ = os.path.dirname(os.path.abspath(__file__))

_amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS',
                'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
                'LEU', 'LYS', 'MET', 'PHE', 'PRO',
                'SER', 'THR', 'TRP', 'TYR', 'VAL']

_nucleotides = ['DA', 'DC', 'DG', 'DT']

class SMOG3SPN2Model(CGModel):
    '''
    A class for SMOG+3SPN2 model.
    '''
    def __init__(self):
        '''
        Initialize. 
        '''
        self.atoms = None
        self.bonded_attr_names = ['protein_bonds', 'protein_angles', 'protein_dihedrals', 'native_pairs', 'exclusions']

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

    
