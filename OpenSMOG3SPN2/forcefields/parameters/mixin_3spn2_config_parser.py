import numpy as np
import pandas as pd
import configparser
import sys
import os

__location__ = os.path.dirname(os.path.abspath(__file__))
dna_3SPN2_conf = f'{__location__}/3SPN2.conf'

_degree_to_rad = np.pi/180
_kcal_to_kj = 4.184
_angstrom_to_nanometer = 0.1

class Mixin3SPN2ConfigParser(object):
    '''
    Define a class for 3SPN2 configuration file parser. 
    Make this an individual class is convenient for other classes to inherit. 
    '''
    def parse_config_file(self, config_file=dna_3SPN2_conf):
        '''
        Parse configuration file. The parameters are loaded as pandas dataframes. 
        Importantly, units in the configuration file 3SPN2.conf are not consistent. 
        This method also converts the units to consistent ones. 
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
        
        self.particle_definition = self.config['Particles']
        self.bond_definition = self.config['Bonds']
        self.angle_definition = self.config['Harmonic Angles']
        self.dihedral_definition = self.config['Dihedrals']
        self.stacking_definition = self.config['Base Stackings']
        self.pair_definition = self.config['Base Pairs']
        self.cross_definition = self.config['Cross Stackings']
        self.protein_dna_particle_definition = self.config['Protein-DNA particles']
        
        # fix units and item names to follow the convention
        # units in 3SPN2.conf are not consistent, and we convert them to consistent units
        # particle definition
        self.particle_definition['epsilon'] *= _kcal_to_kj
        self.particle_definition = self.particle_definition.rename(columns={'radius': 'sigma'})
        self.particle_definition['sigma'] *= _angstrom_to_nanometer
        # bond definition
        self.bond_definition = self.bond_definition.rename(columns={'Kb2': 'k_bond_2', 
                                                                    'Kb3': 'k_bond_3', 
                                                                    'Kb4': 'k_bond_4'})
        # angle definition
        self.angle_definition['k_angle'] = self.angle_definition['epsilon']*2
        self.angle_definition = self.angle_definition.rename(columns={'t0': 'theta0'})
        self.angle_definition['theta0'] *= _degree_to_rad
        # stacking definition
        self.stacking_definition['epsilon'] *= _kcal_to_kj
        self.stacking_definition['sigma'] *= _angstrom_to_nanometer
        self.stacking_definition = self.stacking_definition.rename(columns={'t0': 'theta0'})
        self.stacking_definition['theta0'] *= _degree_to_rad
        self.stacking_definition['alpha'] /= _angstrom_to_nanometer
        # dihedral definition
        self.dihedral_definition['K_dihedral'] *= _kcal_to_kj
        self.dihedral_definition['K_gaussian'] *= _kcal_to_kj
        self.dihedral_definition = self.dihedral_definition.rename(columns={'t0': 'theta0'})
        self.dihedral_definition['theta0'] *= _degree_to_rad
        # base pair definition
        self.pair_definition['torsion'] *= _degree_to_rad
        self.pair_definition['sigma'] *= _angstrom_to_nanometer
        self.pair_definition[['t1', 't2']] *= _degree_to_rad
        self.pair_definition['epsilon'] *= _kcal_to_kj
        self.pair_definition['alpha'] /= _angstrom_to_nanometer
        # cross stacking definition
        self.cross_definition[['t03', 'T0CS_1', 'T0CS_2']] *= _degree_to_rad
        self.cross_definition[['eps_cs1', 'eps_cs2']] *= _kcal_to_kj
        self.cross_definition[['alpha_cs1', 'alpha_cs2']] /= _angstrom_to_nanometer
        self.cross_definition[['Sigma_1', 'Sigma_2']] *= _angstrom_to_nanometer
        # protein-DNA particle definition (parameters for protein-DNA nonbonded interactions)
        self.protein_dna_particle_definition['epsilon'] *= _kcal_to_kj
        self.protein_dna_particle_definition = self.protein_dna_particle_definition.rename(columns={'radius': 'sigma'})
        self.protein_dna_particle_definition[['sigma', 'cutoff']] *= _angstrom_to_nanometer
        


