import numpy as np
import pandas as pd
import sys
import os
import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as unit
import mdtraj

sys.path.append('../..')

from OpenSMOG3SPN2.forcefields.parsers import DNA3SPN2Parser
from OpenSMOG3SPN2.forcefields import SMOG3SPN2Model
from OpenSMOG3SPN2.utils.helper_functions import get_WC_paired_seq

box_a, box_b, box_c = 20.0, 20.0, 20.0
n_dsDNA = 5
dna_type = sys.argv[1]

with open('single_nucl_dna_seq.txt', 'r') as f:
    seq = f.readlines()[0].strip()
full_seq = seq + get_WC_paired_seq(seq)

dna_parser = DNA3SPN2Parser.from_atomistic_pdb('dna.pdb', 'cg_dna.pdb', new_sequence=full_seq, dna_type=dna_type,
                                               temp_name='dna')
top = app.PDBFile('start.pdb').getTopology()
init_coord = app.PDBFile('start.pdb').getPositions()
model = SMOG3SPN2Model(dna_type=dna_type)
for i in range(n_dsDNA):
    model.append_mol(dna_parser)

model.create_system(top, box_a=box_a, box_b=box_b, box_c=box_c)
model.add_dna_bonds(force_group=5)
model.add_dna_angles(force_group=6)
model.add_dna_stackings(force_group=7)
model.add_dna_dihedrals(force_group=8)
model.add_dna_base_pairs(force_group=9)
model.add_dna_cross_stackings(force_group=10)
model.parse_all_exclusions()
model.add_all_vdwl(force_group=11)
model.add_all_elec(force_group=12, salt_conc=100*unit.millimolar)

#model.dna_dihedrals.round(6).to_csv(f'{dna_type}_dna_dihedrals.csv', index=False)

#with open(f'{dna_type}_system.xml', 'w') as f:
#    f.write(mm.XmlSerializer.serialize(model.system))

temperature = 300*unit.kelvin
friction_coeff = 0.01/unit.picosecond
timestep = 10*unit.femtosecond
integrator = mm.LangevinMiddleIntegrator(temperature, friction_coeff, timestep)
model.set_simulation(integrator, platform_name='CUDA', init_coord=init_coord)

traj = mdtraj.load_dcd('traj.dcd', 'start.pdb')
columns = ['dna bond', 'dna angle', 'dna stacking', 'dna dihedral', 'dna base pair', 'dna cross stacking', 
           'all vdwl', 'all elec']
df_energies = pd.DataFrame(columns=columns)
for i in range(traj.xyz.shape[0]):
    model.simulation.context.setPositions(traj.xyz[i])
    row = []
    for j in range(5, 13):
        state = model.simulation.context.getState(getEnergy=True, groups={j})
        energy = state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
        row.append(energy)
    df_energies.loc[len(df_energies.index)] = row
df_energies.round(2).to_csv(f'{dna_type}_energy.csv', index=False)

