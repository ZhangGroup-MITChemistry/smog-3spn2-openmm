import numpy as np
import pandas as pd
import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as unit
import sys
import os

# define some constants based on CODATA
NA = unit.AVOGADRO_CONSTANT_NA # Avogadro constant
kB = unit.BOLTZMANN_CONSTANT_kB  # Boltzmann constant
EC = 1.602176634e-19*unit.coulomb # elementary charge
VEP = 8.8541878128e-12*unit.farad/unit.meter # vacuum electric permittivity

def harmonic_bond_term(df_bonds, use_pbc, force_group=1):
    bonds = mm.HarmonicBondForce()
    for i, row in df_bonds.iterrows():
        a1 = int(row['a1'])
        a2 = int(row['a2'])
        r0 = row['r0']
        k_bond = row['k_bond']
        bonds.addBond(a1, a2, r0, k_bond)
    bonds.setUsesPeriodicBoundaryConditions(use_pbc)
    bonds.setForceGroup(force_group)
    return bonds


def native_pair_gaussian_term(df_native_pairs, use_pbc, force_group=4):
    bonds = mm.CustomBondForce('''energy;
            energy=(-epsilon_G*G+alpha_G*(1-G)/r^12-offset)*step(cutoff-r);
            offset=-epsilon_G*exp(-18)+alpha_G*(1-exp(-18))/cutoff^12;
            cutoff=mu+6*sigma_G;
            G=exp(-(r-mu)^2/(2*sigma_G^2))''')
    bonds.addPerBondParameter('epsilon_G')
    bonds.addPerBondParameter('mu')
    bonds.addPerBondParameter('sigma_G')
    bonds.addPerBondParameter('alpha_G')
    for i, row in df_native_pairs.iterrows():
        a1 = int(row['a1'])
        a2 = int(row['a2'])
        epsilon_G = float(row['epsilon_G'])
        mu = float(row['mu'])
        sigma_G = float(row['sigma_G'])
        alpha_G = float(row['alpha_G'])
        bonds.addBond(a1, a2, [epsilon_G, mu, sigma_G, alpha_G])
    bonds.setUsesPeriodicBoundaryConditions(use_pbc)
    bonds.setForceGroup(force_group)
    return bonds


def native_pair_12_10_term(df_native_pairs, use_pbc, force_group=4):
    '''
    mu is the lowest energy distance for the 12-10 potential.
    '''
    bonds = mm.CustomBondForce('epsilon*(5*(mu/r)^12-6*(mu/r)^10)')
    bonds.addPerBondParameter('epsilon')
    bonds.addPerBondParameter('mu')
    for i, row in df_native_pairs.iterrows():
        a1 = int(row['a1'])
        a2 = int(row['a2'])
        epsilon = row['epsilon']
        mu = row['mu']
        bonds.addBond(a1, a2, [epsilon, mu])
    bonds.setUsesPeriodicBoundaryConditions(use_pbc)
    bonds.setForceGroup(force_group)
    return bonds


def class2_bond_term(df_bonds, use_pbc, force_group=1):
    bonds = mm.CustomBondForce('k_bond_2*(r-r0)^2+k_bond_3*(r-r0)^3+k_bond_4*(r-r0)^4')
    bonds.addPerBondParameter('r0')
    bonds.addPerBondParameter('k_bond_2')
    bonds.addPerBondParameter('k_bond_3')
    bonds.addPerBondParameter('k_bond_4')
    for i, row in df_bonds.iterrows():
        a1 = int(row['a1'])
        a2 = int(row['a2'])
        parameters = row[['r0', 'k_bond_2', 'k_bond_3', 'k_bond_4']].tolist()
        bonds.addBond(a1, a2, parameters)
    bonds.setUsesPeriodicBoundaryConditions(use_pbc)
    bonds.setForceGroup(force_group)
    return bonds



