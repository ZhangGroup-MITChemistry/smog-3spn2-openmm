import numpy as np
import pandas as pd
import sys
import os
import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as unit
import mdtraj

sys.path.append('/home/gridsan/sliu/my-tools/open3spn2')
import open3SPN2

"""
Rerun trajectory with open3SPN2 to compare energy. 
"""

box_a, box_b, box_c = 20.0, 20.0, 20.0

n_dsDNA = 5
dna_type = sys.argv[1]

# directly load CG pdb to ensure P-S-B atom order in each nucleotide
# first load a single dsDNA, then manually combine multiple copies
dsDNA = open3SPN2.DNA.fromCoarsePDB('../test-single-dsDNA/cg_dna.pdb', dna_type=dna_type) 
dsDNA.periodic = True
#dsDNA.dihedrals.round(6).to_csv(f'open3SPN2_{dna_type}_dna_dihedrals.csv', index=False)

multi_dsDNA = open3SPN2.DNA(periodic=True)
multi_dsDNA.parseConfigurationFile()
multi_dsDNA.DNAtype = dna_type
single_dsDNA_atoms = dsDNA.atoms.copy()
multi_dsDNA_atoms = pd.concat([single_dsDNA_atoms]*n_dsDNA, ignore_index=True)
c = 0
old_chainIDs = multi_dsDNA_atoms['chainID'].tolist()
new_chainIDs = []
for i in range(len(multi_dsDNA_atoms.index)):
    if i >= 1:
        if old_chainIDs[i] != old_chainIDs[i - 1]:
            c += 1
    new_chainIDs.append(c)
multi_dsDNA_atoms['chainID'] = new_chainIDs
multi_dsDNA.atoms = multi_dsDNA_atoms

n_atoms_per_dsDNA = len(single_dsDNA_atoms.index)
for i in ['bonds', 'angles', 'stackings', 'dihedrals']:
    a = getattr(dsDNA, i)
    a_list = []
    for j in range(n_dsDNA):
        a_j = a.copy()
        for k in ['aai', 'aaj', 'aak', 'aal']:
            if k in a_j.columns:
                a_j[k] += j*n_atoms_per_dsDNA
        a_list.append(a_j)
    setattr(multi_dsDNA, i, pd.concat(a_list, ignore_index=True))

multi_dsDNA.pdb_file = 'start.pdb'
s = open3SPN2.System(multi_dsDNA, periodicBox=[box_a, box_b, box_c])
s.add3SPN2forces(verbose=True)

#with open(f'open3SPN2_{dna_type}_system.xml', 'w') as f:
#    f.write(mm.XmlSerializer.serialize(s._wrapped_system))

# force group 6-13 for DNA bond, angle, stacking, dihedral, base pair, cross stacking, exclusion, and electrostatics
traj = mdtraj.load_dcd('traj.dcd', 'start.pdb')
columns = ['dna bond', 'dna angle', 'dna stacking', 'dna dihedral', 'dna base pair', 'dna cross stacking', 
           'all vdwl', 'all elec']
s.initializeMD(temperature=300*unit.kelvin, platform_name='CUDA')
simulation = s.simulation
df_energies = pd.DataFrame(columns=columns)
for i in range(traj.xyz.shape[0]):
    simulation.context.setPositions(traj.xyz[i])
    row = []
    for j in range(6, 14):
        state = simulation.context.getState(getEnergy=True, groups={j})
        energy = state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
        row.append(energy)
    df_energies.loc[len(df_energies.index)] = row
df_energies.round(2).to_csv(f'open3SPN2_{dna_type}_energy.csv', index=False)

