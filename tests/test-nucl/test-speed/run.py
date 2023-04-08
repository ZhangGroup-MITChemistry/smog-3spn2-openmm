import numpy as np
import pandas as pd
import sys
import os
import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as unit

sys.path.append('../../..')

from OpenSMOG3SPN2.forcefields.parsers import SMOGParser, DNA3SPN2Parser
from OpenSMOG3SPN2.forcefields import SMOG3SPN2Model
from OpenSMOG3SPN2.utils.helper_functions import get_WC_paired_seq
from OpenSMOG3SPN2.forcefields.rigid import createRigidBodies
from OpenSMOG3SPN2.utils.chromatin_helper_functions import get_chromatin_rigid_bodies

nucl = SMOG3SPN2Model()

# load histone
histone_parser = SMOGParser.from_atomistic_pdb('../pdb-files/histone.pdb', 'histone_CA.pdb', default_parse=False)
histone_parser.parse_mol(get_native_pairs=False)
nucl.append_mol(histone_parser)

# load DNA
with open('dna_seq.txt', 'r') as f:
    seq1 = f.readlines()[0].strip()
seq2 = get_WC_paired_seq(seq1)
target_seq = seq1 + seq2

dna_parser = DNA3SPN2Parser.from_atomistic_pdb('../pdb-files/dna.pdb', new_sequence=target_seq)
nucl.append_mol(dna_parser)
nucl.atoms_to_pdb('cg_nucl.pdb')
nucl.parse_all_exclusions()

top = app.PDBFile('cg_nucl.pdb').getTopology()
init_coord = app.PDBFile('cg_nucl.pdb').getPositions()
rigid_coord = init_coord
nucl.create_system(top, box_a=200, box_b=200, box_c=200)
rigid_bodies = get_chromatin_rigid_bodies(n_nucl=1, nrl=147)
nucl.set_rigid_bodies(rigid_coord, rigid_bodies, keep_unchanged=['protein_exclusions', 'dna_exclusions', 
                                                                 'exclusions'])
#createRigidBodies(nucl.system, rigid_coord, rigid_bodies)
nucl.add_protein_bonds(force_group=1)
nucl.add_protein_angles(force_group=2)
#nucl.add_protein_dihedrals(force_group=3)
#nucl.add_native_pairs(force_group=4)
nucl.add_dna_bonds(force_group=5)
nucl.add_dna_angles(force_group=6)
nucl.add_dna_stackings(force_group=7)
nucl.add_dna_dihedrals(force_group=8)
nucl.add_dna_base_pairs(force_group=9)
nucl.add_dna_cross_stackings(force_group=10)
nucl.add_all_vdwl(force_group=11)
nucl.add_all_elec(force_group=12)

temperature = 300*unit.kelvin
friction_coeff = 0.01/unit.picosecond
timestep = 10*unit.femtosecond
integrator = mm.LangevinMiddleIntegrator(temperature, friction_coeff, timestep)
nucl.set_simulation(integrator, platform_name='CUDA', init_coord=init_coord)
nucl.simulation.minimizeEnergy()
output_interval = 10000
output_dcd = 'output.dcd'
nucl.add_reporters(output_interval, output_dcd)
nucl.simulation.context.setVelocitiesToTemperature(temperature)
nucl.simulation.step(1000000)

    
