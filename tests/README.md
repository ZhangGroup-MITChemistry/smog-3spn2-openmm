# tests

Compare the OpenMM and LAMMPS output energies to validate openmm implementations

Notes:

(1) As LAMMPS read native pairs from SMOG output, here when doing rerun with OpenMM, we use native pairs given by SMOG to keep consistent with LAMMPS settings. 

(2) OpenMM version 3SPN2 has a new set of parameters for building DNA template, while LAMMPS uses the old set of parameters. To recover the old set of parameters with OpenMM, follow the example in `test-tetra-nucl-nrl-172/openmm-rerun-old-geometry/rerun_old_geometry.py`. Such template affects DNA bond and angle potentials. 

(3) Cross stacking interactions are also modified so the cross stacking interactions are slightly different between the outputs of two programs. 

