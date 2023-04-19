# tetra-nucl-nrl-167-conc-0.35mM

To run tetranucleosome simulation in sea of nucleosomes, we need to do following steps: 

(1) Put tetranucleosome and single nucleosomes into a simulation box. 

(2) Fix tetranucleosome (i.e. make the whole tetranucleosome as a rigid body) and relax the coordinates of single nucleosomes with NVT simulation. Here single nuclesome histone core and middle 73 bp core DNA are fixed. 

(3) Release the restriction on tetranucleosome integrity. Just fix histone with middle 73 bp core DNA for tetranucleosome and single nucleosomes and run NVT simulation. 

Practically, we can do the following steps:

(1) Build one openmm system for a given box size and given number of single nucleosomes. In this system, we do not need to add any rigid bodies. 

(2) For selected tetranucleosome configurations, build the initial configuration by inserting tetranucleosome and single nucleosomes. 

(3) Fix tetranucleosome and relax single nucleosomes. 

(4) Run NVT for the whole system. 

