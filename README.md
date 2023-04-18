# smog-3spn2-openmm

An OpenMM implemntation of SMOG and 3SPN2. Note the force fields have been implemented to OpenMM by separate groups. 

The original implementation of OpenMM SMOG (OpenSMOG) by Onuchic group: <https://github.com/smog-server/OpenSMOG>

The original implementation of OpenMM 3SPN2 (open3spn2) by Wolynes group: <https://github.com/cabb99/open3spn2>

## Features

Our implementation has the following advantages: 

(1) All the workflow can be done with Python without requiring SMOG2 outputs, while OpenSMOG requires reading outputs of SMOG2. 

(2) Protein-DNA nonbonded interactions are included, so we can simulate protein-DNA systems. 

(3) It follows OpenABC (a package developed by Zhang group at MIT chemistry) framework, so it is simple to append molecules and many useful functions are included. 

(4) We have optimized 3SPN2 code. Especially, we have optimized the code of cross-stacking interactions. This term is the most time consuming part, and after optimization, for a tetranucleosome system, the speed is about 1.3 times faster than before when running on GPU.

We may combine our OpenSMOG3SPN2 into OpenABC in the future. 

## Notes

Our implementation is intended to be matched with our LAMMPS implementation. Note that in LAMMPS implementation, we keep exclude 1-4 atom pairs from nonbonded interactions even if dihedral potentials are removed (e.g. we removed histone tail dihedrals), while we do not exclude nonbonded interactions between native pair atoms. That is to say, in LAMMPS implementation, the nonbonded exclusion list is composed of all the 1-2, 1-3, and 1-4 pairs. 

