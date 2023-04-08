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
    The class for SMOG+3SPN2 model. 
    To ensure this model works properly, please ensure two neighboring ssDNA chains do not share same chainID. 
    '''
    def __init__(self, dna_type='B_curved', OpenCLPatch=True, default_parse_config=True):
        '''
        Initialize. 
        
        Parameters
        ----------
        dna_type : str
            DNA type. This is related to force field parameters. 
        
        OpenCLPatch : bool
            Whether to use OpenCL patch. 
        
        default_parse_config : bool
            Whether to parse the default 3SPN2 configuration file. 
        
        '''
        self.atoms = None
        self.dna_exclusions = None
        self.exclusions = None
        # note we only set protein_exclusions as bonded attributes, since dna_exclusions and exclusions should be parsed by method after all the molecules are appended
        self.bonded_attr_names = ['protein_bonds', 'protein_angles', 'protein_dihedrals', 'native_pairs', 
                                  'dna_bonds', 'dna_angles', 'dna_stackings', 'dna_dihedrals', 'protein_exclusions']
        self.dna_type = dna_type
        self.OpenCLPatch = OpenCLPatch
        if default_parse_config:
            # load parameters
            self.parse_config_file()
    
    # rewrite append_mol method as we need to check DNA type consistency when appending new DNA
    def append_mol(self, new_mol, verbose=False):
        '''
        The method can append new molecules by concatenating atoms and bonded interaction information saved in dataframes. 
        Please ensure two neighboring chains do not share chainID. 
        
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
        super().append_mol(new_mol, verbose)
    
    def parse_dna_exclusions(self):
        '''
        Parse DNA nonbonded interaction exclusions. 
        Note we parse exclusions here instead of in parser as 3SPN2 involves exclusions between atoms from different chains. 
        The code should be efficient if there are many DNA chains in the system. 
        
        '''
        atoms = self.atoms.copy()
        atoms['index'] = list(range(len(atoms.index)))
        new_chainIDs = []
        c = 0
        for i, row in atoms.iterrows():
            if i >= 1:
                if row['chainID'] != atoms.loc[i - 1, 'chainID']:
                    c += 1
            new_chainIDs.append(c)
        atoms['chainID'] = new_chainIDs
        dna_atoms = atoms[atoms['resname'].isin(_nucleotides)].copy()
        dna_atoms['group'] = dna_atoms['name'].replace(['A', 'T', 'C', 'G'], 'B')
        dna_atoms = dna_atoms.set_index(['chainID', 'resSeq', 'group'])
        
        # set exclusions for atoms from neighboring residues
        dna_exclusions = []
        for i in dna_atoms.index:
            c, r = i[0], i[1]
            a1 = int(dna_atoms.loc[i, 'index'])
            for delta_r in [0, 1]:
                for g in ['P', 'S', 'B']:
                    j = (c, r + delta_r, g)
                    if j in dna_atoms.index:
                        a2 = int(dna_atoms.loc[j, 'index'])
                        if a1 < a2:
                            dna_exclusions.append([a1, a2])

        # set exclusions between W-C base pairs
        if self.OpenCLPatch:
            for b in ['A', 'C']:
                dna_atoms_1 = dna_atoms[dna_atoms['name'] == b]
                dna_atoms_2 = dna_atoms[dna_atoms['name'] == _WC_pair_dict[b]]
                for i in dna_atoms_1['index'].tolist():
                    for j in dna_atoms_2['index'].tolist():
                        a1, a2 = int(i), int(j)
                        if a1 > a2:
                            a1, a2 = a2, a1
                        dna_exclusions.append([a1, a2])
        
        if len(dna_exclusions) >= 1:
            dna_exclusions = pd.DataFrame(np.array(dna_exclusions), columns=['a1', 'a2'])
            dna_exclusions = dna_exclusions.drop_duplicates(ignore_index=True)
            self.dna_exclusions = dna_exclusions.sort_values(by=['a1', 'a2'], ignore_index=True)
        else:
            self.dna_exclusions = pd.DataFrame(columns=['a1', 'a2'])
    
    def parse_all_exclusions(self):
        '''
        Parse all the exclusions (including protein exclusions and DNA exclusions). 
        Run this command before adding nonbonded interactions. 
        
        '''
        self.parse_dna_exclusions()
        if hasattr(self, 'protein_exclusions'):
            if getattr(self, 'protein_exclusions') is None:
                self.protein_exclusions = pd.DataFrame(columns=['a1', 'a2'])
        else:
            self.protein_exclusions = pd.DataFrame(columns=['a1', 'a2'])
        exclusions = pd.concat([self.protein_exclusions, self.dna_exclusions], ignore_index=True)
        self.exclusions = exclusions.sort_values(by=['a1', 'a2'], ignore_index=True)
    
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
        
        Please ensure two neighboring chains in self.atoms do not share the same chainID, so that different chains can be distinguished properly. 
        
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
        atoms.index = list(range(len(atoms.index)))
        new_chainIDs = []
        c = 0
        for i, row in atoms.iterrows():
            if i >= 1:
                if row['chainID'] != atoms.loc[i - 1, 'chainID']:
                    c += 1
            new_chainIDs.append(c)
        atoms['chainID'] = new_chainIDs
        atoms['index'] = list(range(len(atoms.index)))
        # originally index is set as tuples of chainID, resSeq, and name
        # we use multiindex instead, which is the more canonical method
        atoms = atoms.set_index(['chainID', 'resSeq', 'name']) # use multiple columns as index
        dna_atoms = atoms[atoms['resname'].isin(_nucleotides)].copy()
        # find W-C pairs
        # note as we have two types of W-C pairs (A-T and C-G), we need two separate forces
        for i, row in pair_definition.iterrows():
            parameters = [row['torsion'], row['sigma'], row['t1'], row['t2'], row['rang'], 
                          row['epsilon'], row['alpha']]
            base1, base2 = row['Base1'], row['Base2']
            donors1 = dna_atoms.loc[pd.IndexSlice[:, :, [base1]]].copy()
            acceptors1 = dna_atoms.loc[pd.IndexSlice[:, :, [base2]]].copy()
            donors2 = dna_atoms.loc[[(x[0], x[1], 'S') for x in donors1.index]].copy()
            acceptors2 = dna_atoms.loc[[(x[0], x[1], 'S') for x in acceptors1.index]].copy()
            force_i = functional_terms.dna_3spn2_base_pair_term(self.use_pbc, cutoff, force_group)
            # add donors and acceptors
            for a1, a2 in zip(donors1['index'].tolist(), donors2['index'].tolist()):
                force_i.addDonor(a1, a2, -1, parameters)
            for a1, a2 in zip(acceptors1['index'].tolist(), acceptors2['index'].tolist()):
                force_i.addAcceptor(a1, a2, -1)
            # set exclusions
            # for hbond forces, exclusions are added based on donor id and acceptor id, instead of atom id
            donors1['donor_id'] = list(range(len(donors1.index)))
            acceptors1['acceptor_id'] = list(range(len(acceptors1.index)))
            for j in donors1.index:
                c, r = j[0], j[1]
                for delta_r in [-2, -1, 1, 2]:
                    k = (c, r + delta_r, base2)
                    if k in acceptors1.index:
                        force_i.addExclusion(int(donors1.loc[j, 'donor_id']), int(acceptors1.loc[k, 'acceptor_id']))
            self.system.addForce(force_i)
    
    def add_dna_cross_stackings(self, cutoff=1.8, force_group=10):
        '''
        Add DNA cross stacking potentials. 
        
        Please ensure two neighboring chains in self.atoms do not share the same chainID, so that different chains can be distinguished properly. 
        
        The method should be efficient if OpenCLPatch is True. 
        
        Parameters
        ----------
        cutoff : float or int
            Cutoff distance. 
        
        force_group : int
            Force group.
        
        '''
        print('Add DNA cross stackings.')
        cross_definition = self.cross_definition[self.cross_definition['DNA'] == self.dna_type].copy()
        cross_definition = cross_definition.set_index(['Base_d1', 'Base_a1', 'Base_a3'])
        atoms = self.atoms.copy()
        # reset chainID to unique numbers so forces can be set properly
        atoms.index = list(range(len(atoms.index)))
        new_chainIDs = []
        c = 0
        for i, row in atoms.iterrows():
            if i >= 1:
                if row['chainID'] != atoms.loc[i - 1, 'chainID']:
                    c += 1
            new_chainIDs.append(c)
        atoms['chainID'] = new_chainIDs
        atoms['index'] = list(range(len(atoms.index)))
        dna_atoms = atoms[atoms['resname'].isin(_nucleotides)].copy()
        dna_atoms['group'] = dna_atoms['name'].replace(['A', 'T', 'C', 'G'], 'B')
        dna_atoms = dna_atoms.set_index(['chainID', 'resSeq', 'group'])
        bases = dna_atoms.loc[pd.IndexSlice[:, :, ['B']]].copy()
        # define forces
        dict_cross_stackings = {}
        for b in ['A', 'T', 'G', 'C']:
            force1 = functional_terms.dna_3spn2_cross_stacking_term(self.use_pbc, cutoff, force_group)
            force2 = functional_terms.dna_3spn2_cross_stacking_term(self.use_pbc, cutoff, force_group)
            dict_cross_stackings.update({b: (force1, force2)})
        donor_dict = {i: [] for i in ['A', 'T', 'G', 'C']}
        acceptor_dict = {i: [] for i in ['A', 'T', 'G', 'C']}
        for i in bases.index:
            a1, a1n = int(bases.loc[i, 'index']), bases.loc[i, 'name']
            a2 = int(dna_atoms.loc[(i[0], i[1], 'S'), 'index'])
            j = (i[0], i[1] + 1, 'B')
            if j in bases.index:
                a3, a3n = int(bases.loc[j, 'index']), bases.loc[j, 'name']
                force1, force2 = dict_cross_stackings[a1n]
                columns = ['t03', 'T0CS_2', 'rng_cs2', 'rng_bp', 'eps_cs2', 'alpha_cs2', 'Sigma_2']
                p = cross_definition.loc[(_WC_pair_dict[a1n], a1n, a3n), columns].tolist()
                force1.addDonor(a1, a2, a3)
                force2.addAcceptor(a1, a2, a3, p)
                donor_dict[a1n].append(a1) # add force1 donor a1
            k = (i[0], i[1] - 1, 'B')
            if k in bases.index:
                a3, a3n = int(bases.loc[k, 'index']), bases.loc[k, 'name']
                force1, force2 = dict_cross_stackings[_WC_pair_dict[a1n]]
                columns = ['t03', 'T0CS_1', 'rng_cs1', 'rng_bp', 'eps_cs1', 'alpha_cs1', 'Sigma_1']
                p = cross_definition.loc[(_WC_pair_dict[a1n], a1n, a3n), columns].tolist()
                force1.addAcceptor(a1, a2, a3, p)
                force2.addDonor(a1, a2, a3)
                acceptor_dict[_WC_pair_dict[a1n]].append(a1) # add force1 acceptor a1
        # set exclusions
        for b in ['A', 'T', 'C', 'G']:
            force1, force2 = dict_cross_stackings[b]
            donors = dna_atoms[dna_atoms['index'].isin(donor_dict[b])].copy()
            acceptors = dna_atoms[dna_atoms['index'].isin(acceptor_dict[b])].copy()
            donors['donor_id'] = list(range(len(donors.index)))
            acceptors['acceptor_id'] = list(range(len(acceptors.index)))
            if self.OpenCLPatch:
                # use more efficient method instead of looping over all the donors and acceptors
                max_delta_resSeq = 2
                for i in donors.index:
                    donor_id = int(donors.loc[i, 'donor_id'])
                    for delta_resSeq in range(-1*max_delta_resSeq, max_delta_resSeq + 1):
                        j = (i[0], i[1] + delta_resSeq, 'B')
                        if j in acceptors.index:
                            acceptor_id = int(acceptors.loc[j, 'acceptor_id'])
                            force1.addExclusion(donor_id, acceptor_id)
                            force2.addExclusion(acceptor_id, donor_id)
            else:
                # loop over all the donors and acceptors
                # this is inefficient, but usually we set self.OpenCLPatch as True and we do not use this method
                max_delta_resSeq = 3
                for i in donors.index:
                    for j in acceptors.index:
                        a1 = int(donors.loc[i, 'index'])
                        a2 = int(acceptors.loc[j, 'index'])
                        donor_id = int(donors.loc[i, 'donor_id'])
                        acceptor_id = int(acceptor_id.loc[j, 'acceptor_id'])
                        if (i[0] == j[0]) and (abs(i[1] - j[1]) <= max_delta_resSeq) or (a1 > a2):
                            # Question: is this correct? It looks weird that as long as a1 > a2 there is an exclusion. 
                            force1.addExclusion(donor_id, acceptor_id)
                            force2.addExclusion(acceptor_id, donor_id)
            self.system.addForce(force1)
            self.system.addForce(force2)
    
    def add_all_vdwl(self, param_PP_MJ_path=f'{__location__}/parameters/pp_MJ.csv', cutoff_PD=1.425*unit.nanometer, 
                     force_group=11):
        '''
        CG atom type 0-19 for amino acids.
        CG atom type 20-25 for DNA atoms. 
        
        Parameters
        ----------
        param_PP_MJ : str
            Protein-protein MJ potential parameter file path. 
        
        '''
        print('Add all the nonbonded contact interactions.')
        param_PP_MJ = pd.read_csv(param_PP_MJ_path)
        force = functional_terms.all_smog_MJ_3spn2_term(self, param_PP_MJ, cutoff_PD, force_group)
        self.system.addForce(force)
    
    def add_all_elec(self, salt_concentration=150*unit.millimolar, temperature=300*unit.kelvin, 
                     elec_DD_charge_scale=0.6, cutoff_DD=5*unit.nanometer, cutoff_PP_PD=3.141504539*unit.nanometer, 
                     dielectric_PP_PD=78, force_group=12):
        print('Add all the electrostatic interactions.')
        force = functional_terms.all_smog_3spn2_elec_term(self, salt_concentration, temperature, 
                                                          elec_DD_charge_scale, cutoff_DD, cutoff_PP_PD, 
                                                          dielectric_PP_PD, force_group)
        self.system.addForce(force)
        

    
