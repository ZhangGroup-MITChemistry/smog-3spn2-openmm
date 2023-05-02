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

box_a, box_b, box_c = 25.0, 25.0, 25.0

n_dsDNA = 1
dna_type = sys.argv[1]

# directly load CG pdb to ensure P-S-B atom order in each nucleotide
dsDNA = open3SPN2.DNA.fromCoarsePDB('cg_dna.pdb', dna_type=dna_type) 
dsDNA.periodic = True
#dsDNA.dihedrals.round(6).to_csv(f'open3SPN2_{dna_type}_dna_dihedrals.csv', index=False)

s = open3SPN2.System(dsDNA, periodicBox=[box_a, box_b, box_c])
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
df_energies.round(4).to_csv(f'open3SPN2_{dna_type}_energy.csv', index=False)

