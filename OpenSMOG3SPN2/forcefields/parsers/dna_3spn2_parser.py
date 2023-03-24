import numpy as np
import pandas as pd
import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as unit
from OpenSMOG3SPN2.utils import helper_functions
import configparser
import sys
import os

__location__ = os.path.dirname(os.path.abspath(__file__))

dna_3SPN2_conf = f'{__location__}/../parameters/3SPN2.conf'

_masses = {"H": 1.00794, "C": 12.0107, "N": 14.0067, "O": 15.9994, "P": 30.973762} 

_CG_map = {"O5\'": 'P', "C5\'": 'S', "C4\'": 'S', "O4\'": 'S', "C3\'": 'S', "O3\'": 'P',
           "C2\'": 'S', "C1\'": 'S', "O5*": 'P', "C5*": 'S', "C4*": 'S', "O4*": 'S',
           "C3*": 'S', "O3*": 'P', "C2*": 'S', "C1*": 'S', "N1": 'B', "C2": 'B', "O2": 'B',
           "N2": 'B', "N3": 'B', "C4": 'B', "N4": 'B', "C5": 'B', "C6": 'B', "N9": 'B',
           "C8": 'B', "O6": 'B', "N7": 'B', "N6": 'B', "O4": 'B', "C7": 'B', "P": 'P',
           "OP1": 'P', "OP2": 'P', "O1P": 'P', "O2P": 'P', "OP3": 'P', "HO5'": 'P',
           "H5'": 'S', "H5''": 'S', "H4'": 'S', "H3'": 'S', "H2'": 'S', "H2''": 'S',
           "H1'": 'S', "H8": 'B', "H61": 'B', "H62": 'B', 'H2': 'B', 'H1': 'B', 'H21': 'B',
           'H22': 'B', 'H3': 'B', 'H71': 'B', 'H72': 'B', 'H73': 'B', 'H6': 'B', 'H41': 'B',
           'H42': 'B', 'H5': 'B', "HO3'": 'P'}

'''
Open3SPN2 was originally developed by Carlos Bueno. 
Most code is adapted from the original Open3SPN2. 
'''


class DNA3SPN2Parser():
    '''
    DNA 3SPN2 parser. 
    '''
    def __init__(self, cg_pdb, seq=None, default_parse=True):
        '''
        Initialize DNA with 3SPN2 model. 
        '''
        pass # to be continued
        # update sequence
        # parse molecule
    
    
    def parse_config_file(self, config_file=dna_3SPN2_conf):
        '''
        Parse configuration file. The parameters are loaded as pandas dataframes. 
        '''
        
        def parse_row(row):
            '''
            Parse one row from the configuration. 
            '''
            values = row.split('#')[0].split() # remove comments
            for i in range(len(values)):
                values[i] = values[i].strip()
                try: 
                    x = int(values[i])
                except ValueError:
                    try:
                        x = float(values[i])
                    except ValueError:
                        x = values[i]
                values[i] = x
            return values
        
        config = configparser.ConfigParser()
        print(config_file)
        config.read(config_file, encoding='utf-8')
        self.config = {}
        for i in config.sections():
            data = []
            for j in config[i]:
                if j == 'name':
                    columns = parse_row(config[i][j])
                elif len(j) > 3 and j[:3] == 'row':
                    data += [parse_row(config[i][j])]
            self.config[i] = pd.DataFrame(data, columns=columns)
    
    
    def parse_mol(self):
        '''
        Parse molecule. 
        '''
        pass
    
    
    @staticmethod
    def aa_to_cg(df_pdb, lammps_format=True):
        '''
        Convert DNA all-atom structure to CG structure. 
        Both input and output structures are saved as pandas dataframe. 
        
        Parameters
        ----------
        df_pdb : pd.DataFrame
            Input all-atom structure. 
        
        lammps_format : bool
            Ensure CG atom order in each nucleotide is P-S-B, while the first nucleotide does not have P. 
        
        Returns
        -------
        
        
        '''
        columns = ['recname', 'serial', 'name', 'altLoc',
                    'resname', 'chainID', 'resSeq', 'iCode',
                    'x', 'y', 'z', 'occupancy', 'tempFactor',
                    'element', 'charge', 'type']
        temp = df_pdb.copy()
        temp = temp[temp['resname'].isin(['DA', 'DT', 'DG', 'DC'])] # select DNA
        temp['group'] = temp['name'].replace(_CG_map) # assign each atom to phosphate, sugar, or base
        temp = temp[temp['group'].isin(['P', 'S', 'B'])]
        
        # Move the O3' to the next residue
        # also remove the O3' atom at the final residue of the chain
        for c in temp['chainID'].unique():
            sel = temp.loc[(temp['name'] == "O3\'") & (temp['chainID'] == c), "resSeq"]
            temp.loc[(temp['name'] == "O3\'") & (temp['chainID'] == c), "resSeq"] = list(sel)[1:] + [-1]
            sel = temp.loc[(temp['name'] == "O3\'") & (temp['chainID'] == c), "resname"]
            temp.loc[(temp['name'] == "O3\'") & (temp['chainID'] == c), "resname"] = list(sel)[1:] + ["remove"]
        #temp = temp[temp['resSeq'] > 0]
        temp = temp[temp['resname'] != 'remove']

        # Calculate center of mass
        # perform multiplication with numpy array
        temp['element'] = temp['element'].str.strip() # remove white spaces on both ends
        temp['mass'] = temp.element.replace(_masses).astype(float)
        coord = temp[['x', 'y', 'z']].to_numpy()
        weight = temp['mass'].to_numpy()
        weighted_coord = (coord.T*weight).T
        temp['mass'] = weighted_coord
        temp = temp[temp['element'] != 'H']  # Exclude hydrogens
        Coarse = temp.groupby(['chainID', 'resSeq', 'resname', 'group']).sum(numeric_only=True).reset_index()
        weighted_coord = Coarse[['x', 'y', 'z']].to_numpy()
        weight = Coarse['mass'].to_numpy()
        Coarse[['x', 'y', 'z']] = (weighted_coord.T/weight).T

        # Set pdb columns
        Coarse['recname'] = 'ATOM'
        Coarse['name'] = Coarse['group']
        Coarse['altLoc'] = ''
        Coarse['iCode'] = ''
        Coarse['occupancy'] = ''
        Coarse['tempFactor'] = ''
        Coarse['charge'] = ''
        # Change name of base to real base
        mask = (Coarse['name'] == 'B')
        Coarse.loc[mask, 'name'] = Coarse[mask].resname.str[-1]  # takes last letter from the residue name
        Coarse['type'] = Coarse['name']
        # Set element (depends on base)
        Coarse['element'] = Coarse['name'].replace({'P': 'P', 'S': 'H', 'A': 'N', 'T': 'S', 'G': 'C', 'C': 'O'})
        # Remove P from the beggining
        drop_list = []
        for chain in Coarse.chainID.unique():
            sel = Coarse[Coarse.chainID == chain]
            drop_list += list(sel[(sel.resSeq == sel.resSeq.min()) & sel['name'].isin(['P'])].index)
        Coarse = Coarse.drop(drop_list)
        
        if lammps_format:
            # rearrange CG atom order
            Coarse['chainID_resSeq'] = Coarse['chainID'].astype(str) + '_' + Coarse['resSeq'].astype(str)
            Coarse.index = Coarse['chainID_resSeq'].astype(str) + '_' + Coarse['group'].astype(str)
            chainID_resSeq_no_duplicates = Coarse['chainID_resSeq'].drop_duplicates(keep='first').tolist()
            Coarse_new = pd.DataFrame(columns=list(Coarse.columns))
            for each in chainID_resSeq_no_duplicates:
                each_P = f'{each}_P'
                if each_P in Coarse.index:
                    Coarse_new.loc[len(Coarse_new.index)] = Coarse.loc[each_P]
                Coarse_new.loc[len(Coarse_new.index)] = Coarse.loc[f'{each}_S']
                Coarse_new.loc[len(Coarse_new.index)] = Coarse.loc[f'{each}_B']
            # check if the number of DNA atoms are the same
            assert len(Coarse_new.index) == len(Coarse.index)
            Coarse = Coarse_new.copy()
        
        # Renumber
        Coarse.index = list(range(len(Coarse.index)))
        Coarse['serial'] = Coarse.index
        return Coarse[columns]


