import numpy as np
import pandas as pd
import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as unit
import math
import sys
import os

_amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS',
                'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
                'LEU', 'LYS', 'MET', 'PHE', 'PRO',
                'SER', 'THR', 'TRP', 'TYR', 'VAL']

_nucleotides = ['DA', 'DT', 'DC', 'DG']

_dna_atom_names = ['P', 'S', 'A', 'T', 'C', 'G']

# define some constants based on CODATA
NA = unit.AVOGADRO_CONSTANT_NA # Avogadro constant
kB = unit.BOLTZMANN_CONSTANT_kB  # Boltzmann constant
EC = 1.602176634e-19*unit.coulomb # elementary charge
VEP = 8.8541878128e-12*unit.farad/unit.meter # vacuum electric permittivity

def all_smog_MJ_3spn2_term(mol, param_PP_MJ, cutoff_PD=1.425*unit.nanometer, force_group=11):
    '''
    Combine all the SMOG (MJ potential for protein-protein nonbonded pairs) and 3SPN2 nonbonded interactions into one force.
    CG atom type 0-19 for amino acids.
    CG atom type 20-25 for DNA atoms. 
    Many parameters are saved as attributes of mol. 
        
    '''
    vdwl = mm.CustomNonbondedForce('''energy;
           energy=4*epsilon*((sigma/r)^12-(sigma/r)^6-offset)*step(cutoff-r);
           offset=(sigma/cutoff)^12-(sigma/cutoff)^6;
           epsilon=epsilon_map(atom_type1, atom_type2);
           sigma=sigma_map(atom_type1, atom_type2);
           cutoff=cutoff_map(atom_type1, atom_type2)''')
    n_atom_types = len(_amino_acids) + len(_dna_atom_names)
    epsilon_map = np.zeros((n_atom_types, n_atom_types))
    sigma_map = np.zeros((n_atom_types, n_atom_types))
    cutoff_map = np.zeros((n_atom_types, n_atom_types))
    # add protein-protein interactions
    for _, row in param_PP_MJ.iterrows():
        atom_type1, atom_type2 = row['atom_type1'], row['atom_type2']
        i = _amino_acids.index(atom_type1)
        j = _amino_acids.index(atom_type2)
        epsilon_map[i, j] = row['epsilon (kj/mol)']
        epsilon_map[j, i] = epsilon_map[i, j]
        sigma_map[i, j] = row['sigma (nm)']
        sigma_map[j, i] = sigma_map[i, j]
        cutoff_map[i, j] = row['cutoff_LJ (nm)']
        cutoff_map[j, i] = cutoff_map[i, j]
    # add DNA-DNA interactions
    param_DD = mol.particle_definition[mol.particle_definition['DNA'] == mol.dna_type].copy()
    param_DD.index = param_DD['name'] # rearrange to make sure the row order is based on dna_atom_names
    param_DD = param_DD.loc[_dna_atom_names]
    param_DD.index = list(range(len(param_DD.index)))
    for i1 in range(len(_dna_atom_names)):
        for j1 in range(i1, len(_dna_atom_names)):
            i = i1 + len(_amino_acids)
            j = j1 + len(_amino_acids)
            epsilon_i = param_DD.loc[i1, 'epsilon']
            epsilon_j = param_DD.loc[j1, 'epsilon']
            epsilon_map[i, j] = (epsilon_i*epsilon_j)**0.5
            epsilon_map[j, i] = epsilon_map[i, j]
            sigma_i = param_DD.loc[i1, 'sigma']
            sigma_j = param_DD.loc[j1, 'sigma']
            sigma_map[i, j] = 0.5*(sigma_i + sigma_j)*(2**(-1/6))
            sigma_map[j, i] = sigma_map[i, j]
            cutoff_map[i, j] = 0.5*(sigma_i + sigma_j)
            cutoff_map[j, i] = cutoff_map[i, j]
    # add protein-DNA interactions
    all_param_PD = mol.protein_dna_particle_definition
    param_dna_PD = all_param_PD[(all_param_PD['molecule'] == 'DNA') & (all_param_PD['DNA'] == mol.dna_type)].copy()
    param_dna_PD.index = param_dna_PD['name']
    param_dna_PD = param_dna_PD.loc[_dna_atom_names].copy()
    param_dna_PD.index = list(range(len(param_dna_PD.index)))
    param_protein_PD =  all_param_PD[(all_param_PD['molecule'] == 'Protein')].copy()
    param_protein_PD.index = param_protein_PD['name']
    param_protein_PD = param_protein_PD.loc[['CA']].copy() # protein only has CA type CG atom
    param_protein_PD.index = list(range(len(param_protein_PD.index)))
    param_PD = pd.concat([param_dna_PD, param_protein_PD], ignore_index=True)
    for i1 in range(len(_dna_atom_names)):
        i = i1 + len(_amino_acids)
        epsilon_i = param_PD.loc[i1, 'epsilon']
        epsilon_j = param_PD.loc[len(_dna_atom_names), 'epsilon']
        epsilon_map[i, :len(_amino_acids)] = (epsilon_i*epsilon_j)**0.5
        epsilon_map[:len(_amino_acids), i] = epsilon_map[i, :len(_amino_acids)]
        sigma_i = param_PD.loc[i1, 'sigma']
        sigma_j = param_PD.loc[len(_dna_atom_names), 'sigma']
        sigma_map[i, :len(_amino_acids)] = 0.5*(sigma_i + sigma_j)
        sigma_map[:len(_amino_acids), i] = sigma_map[i, :len(_amino_acids)]
        cutoff_map[i, :len(_amino_acids)] = cutoff_PD.value_in_unit(unit.nanometer)
        cutoff_map[:len(_amino_acids), i] = cutoff_map[i, :len(_amino_acids)]
    max_cutoff = np.amax(cutoff_map)
    epsilon_map = epsilon_map.ravel().tolist()
    sigma_map = sigma_map.ravel().tolist()
    cutoff_map = cutoff_map.ravel().tolist()
    vdwl.addTabulatedFunction('epsilon_map', mm.Discrete2DFunction(n_atom_types, n_atom_types, epsilon_map))
    vdwl.addTabulatedFunction('sigma_map', mm.Discrete2DFunction(n_atom_types, n_atom_types, sigma_map))
    vdwl.addTabulatedFunction('cutoff_map', mm.Discrete2DFunction(n_atom_types, n_atom_types, cutoff_map))
    vdwl.addPerParticleParameter('atom_type')
    # add atom type
    for i, row in mol.atoms.iterrows():
        resname_i = row['resname']
        name_i = row['name']
        if (resname_i in _amino_acids) and (name_i == 'CA'):
            vdwl.addParticle([_amino_acids.index(resname_i)])
        elif (resname_i in _nucleotides) and (name_i in _dna_atom_names):
            vdwl.addParticle([len(_amino_acids) + _dna_atom_names.index(name_i)])
        else:
            sys.exit(f'Cannot recognize atom with resname {resname_i} and name {name_i}.')
    # add exclusions
    for i, row in mol.exclusions.iterrows():
        vdwl.addExclusion(int(row['a1']), int(row['a2']))
    # set PBC, cutoff, and force group
    if mol.use_pbc:
        vdwl.setNonbondedMethod(vdwl.CutoffPeriodic)
    else:
        vdwl.setNonbondedMethod(vdwl.CutoffNonPeriodic)
    vdwl.setCutoffDistance(max_cutoff)
    vdwl.setForceGroup(force_group)
    return vdwl


def all_smog_3spn2_elec_term(mol, salt_conc=150*unit.millimolar, temperature=300*unit.kelvin, 
                             elec_DD_charge_scale=0.6, cutoff_DD=5*unit.nanometer, 
                             cutoff_PP_PD=3.141504539*unit.nanometer, dielectric_PP_PD=78, force_group=12):
    '''
    Combine all the SMOG and 3SPN2 electrostatic interactions into one force. 
    
    CG atom types: 
    Type 0 for zero-charge CG atoms. 
    Type 1 for ARG and LYS CA atoms. 
    Type 2 for ASP and GLU CA atoms. 
    Type 3 for phosphate CG atoms. 

    '''
    C = salt_conc.value_in_unit(unit.molar)
    T = temperature.value_in_unit(unit.kelvin)
    print(f'For electrostatic interactions, set monovalent salt concentration as {1000*C} mM.')
    print(f'For electrostatic interactions, set temperature as {T} K.')
    e = 249.4 - 0.788*T + 7.2E-4*T**2
    a = 1 - 0.2551*C + 5.151E-2*C**2 - 6.889E-3*C**3
    dielectric_DD = e*a
    print(f'DNA-DNA dielectric constant is {dielectric_DD}')
    print(f'Protein-protein and protein-DNA dielectric constant is {dielectric_PP_PD}.')
    elec = mm.CustomNonbondedForce('''energy;
           energy=alpha*exp(-r/ldby)*step(cutoff-r)/r;
           alpha=alpha_map(cg_atom_type1, cg_atom_type2);
           ldby=ldby_map(cg_atom_type1, cg_atom_type2);
           cutoff=cutoff_map(cg_atom_type1, cg_atom_type2)''')
    n_atom_types = 4
    charge_list = [0, 1, -1, -1]
    # use Discrete2DFunction to define mappings for alpha, sigma, and cutoff
    alpha_map = np.zeros((n_atom_types, n_atom_types))
    ldby_map = np.zeros((n_atom_types, n_atom_types))
    cutoff_map = np.zeros((n_atom_types, n_atom_types))
    # set mappings
    for i in range(n_atom_types):
        for j in range(i, n_atom_types):
            q_i = charge_list[i]
            q_j = charge_list[j]
            if (i == 3) and (j == 3):
                # phosphate-phosphate electrostatic interactions
                q_i *= elec_DD_charge_scale
                q_j *= elec_DD_charge_scale
                cutoff_ij = cutoff_DD
                denominator = 4*np.pi*VEP*dielectric_DD/(NA*(EC**2))
                denominator = denominator.value_in_unit(unit.kilojoule_per_mole**-1*unit.nanometer**-1)
                ldby_ij = (dielectric_DD*VEP*kB*temperature/(2.0*NA*(EC**2)*salt_conc))**0.5
            else:
                cutoff_ij = cutoff_PP_PD
                denominator = 4*np.pi*VEP*dielectric_PP_PD/(NA*(EC**2))
                denominator = denominator.value_in_unit(unit.kilojoule_per_mole**-1*unit.nanometer**-1)
                ldby_ij = (dielectric_PP_PD*VEP*kB*temperature/(2.0*NA*(EC**2)*salt_conc))**0.5
            alpha_map[i, j] = q_i*q_j/denominator
            alpha_map[j, i] = alpha_map[i, j]
            cutoff_map[i, j] = cutoff_ij.value_in_unit(unit.nanometer)
            cutoff_map[j, i] = cutoff_map[i, j]
            ldby_map[i, j] = ldby_ij.value_in_unit(unit.nanometer)
            ldby_map[j, i] = ldby_map[i, j]
    max_cutoff = np.amax(cutoff_map)
    alpha_map = alpha_map.ravel().tolist()
    ldby_map = ldby_map.ravel().tolist()
    cutoff_map = cutoff_map.ravel().tolist()
    elec.addTabulatedFunction('alpha_map', mm.Discrete2DFunction(n_atom_types, n_atom_types, alpha_map))
    elec.addTabulatedFunction('ldby_map', mm.Discrete2DFunction(n_atom_types, n_atom_types, ldby_map))
    elec.addTabulatedFunction('cutoff_map', mm.Discrete2DFunction(n_atom_types, n_atom_types, cutoff_map))
    elec.addPerParticleParameter('cg_atom_type')
    # add atom type
    for i, row in mol.atoms.iterrows():
        resname_i = row['resname']
        name_i = row['name']
        if (resname_i in ['ARG', 'LYS']) and (name_i == 'CA'):
            elec.addParticle([1])
        elif (resname_i in ['ASP', 'GLU']) and (name_i == 'CA'):
            elec.addParticle([2])
        elif (resname_i in _nucleotides) and (name_i == 'P'):
            elec.addParticle([3])
        else:
            elec.addParticle([0])
    # add exclusions
    for i, row in mol.exclusions.iterrows():
        elec.addExclusion(int(row['a1']), int(row['a2']))
    # set PBC, cutoff, and force group
    if mol.use_pbc:
        elec.setNonbondedMethod(elec.CutoffPeriodic)
    else:
        elec.setNonbondedMethod(elec.CutoffNonPeriodic)
    elec.setCutoffDistance(max_cutoff)
    elec.setForceGroup(force_group)
    return elec



