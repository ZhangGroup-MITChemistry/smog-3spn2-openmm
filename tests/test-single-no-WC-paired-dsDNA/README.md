# test-single-no-WC-paired-dsDNA

Rerun the trajectory with not W-C paired sequence and B_curved DNA model. This will enforce OpenSMOG3SPN2 to use the same DNA template as open3SPN2. In principle this should make open3SPN2 and OpenSMOG3SPN2 to give identical potential energies for DNA bond, angle, and dihedrals. 

We just use the trajectory in test-single-dsDNA. Both reruns are using GPU with "mixed" precision. 

