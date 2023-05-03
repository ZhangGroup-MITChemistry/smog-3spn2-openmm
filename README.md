# smog-3spn2-openmm

An OpenMM implemntation of SMOG and 3SPN2. Note the force fields have been implemented to OpenMM by separate groups. 

The original implementation of OpenMM SMOG (OpenSMOG) by Onuchic group: <https://github.com/smog-server/OpenSMOG>

The original implementation of OpenMM 3SPN2 (open3spn2) by Wolynes group: <https://github.com/cabb99/open3spn2>

## Features

Our implementation has the following features: 

(1) All the workflow can be done with Python without requiring SMOG2 outputs, while OpenSMOG requires reading outputs of SMOG2. 

(2) Protein-DNA nonbonded interactions are included, so we can simulate protein-DNA systems. 

(3) It follows OpenABC (a package developed by Zhang group at MIT chemistry) framework, so it is simple to append molecules and many useful functions are included. This makes simulating multiple DNA chains much easier. 

(4) We have optimized 3SPN2 code performance. Especially, we have optimized the code of cross-stacking interactions. This term is the most time consuming part, and after optimization, for a tetranucleosome system, the speed is about 1.3 times faster than before when running on GPU.

(5) The original open3SPN2 code has two modes: OpenCLPatch as True or False. Our implementation here is equivalent to open3SPN2 OpenCLPatch = True mode. 

We may combine our OpenSMOG3SPN2 into OpenABC in the future. 

## Manual

Manual is at: <https://zhanggroup-mitchemistry.github.io/smog-3spn2-openmm/>

## Notes

Our implementation is intended to be matched with our LAMMPS implementation. Note that in LAMMPS implementation, we keep exclude 1-4 atom pairs from nonbonded interactions even if dihedral potentials are removed (e.g. we removed histone tail dihedrals), while we do not exclude nonbonded interactions between native pair atoms. That is to say, in LAMMPS implementation, the nonbonded exclusion list is composed of all the 1-2, 1-3, and 1-4 pairs. 

We have only compared 3SPN2 "B_curved" DNA type energy with LAMMPS, while A and B-DNA are not compared with LAMMPS. Meanwhile, we only compare "OpenCLPatch" as True mode. A-DNA and B-DNA results are compared with open3SPN2 results. 

## Citations

If you found this useful, please consider citing the following papers: 

*Open3SPN2*: Lu, Wei, et al. "OpenAWSEM with Open3SPN2: A fast, flexible, and accessible framework for large-scale coarse-grained biomolecular simulations." PLoS computational biology 17.2 (2021): e1008308.

*OpenSMOG*: de Oliveira Jr, Antonio B., et al. "SMOG 2 and OpenSMOG: Extending the limits of structure‐based models." Protein Science 31.1 (2022): 158-172.

*OpenABC*: Liu, Shuming, et al. "OpenABC Enables Flexible, Simplified, and Efficient GPU Accelerated Simulations of Biomolecular Condensates." doi: https://doi.org/10.1101/2023.04.19.537533

