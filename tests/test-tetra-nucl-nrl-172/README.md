Compare OpenMM and LAMMPS output energies for tetranucleosome. 

pairs.dat includes native pairs found by SMOG. The columns are: aai, aaj, type, epsilon_G, mu, sigma_G, and alpha_G. 

ca_pair_list.txt includes native pairs used by LAMMPS. 

Note the native pairs for each histone found by SMOG are not identical. The good news is we use rigid bodies so native pairs do not affect simulations. 

Important: in LAMMPS, even though we remove dihedrals involving histone tail atoms, and we keep native pairs within histone cores, nonbonded interactions between 1-4 atoms involving histone tail atoms are excluded, and nonbonded interactions between native pair atoms are not excluded!

