# test-single-DNA

We test our code for single DNA to make sure the modified code gives identical energy to the original open3SPN2. 

First use run_md.py to produce a trajectory for multiple nucleosomal dsDNA (histones are not included), then use rerun.py to rerun the trajectory with OpenSMOG3SPN2 and open3SPN2 to compare. We need to input dna_type, and here in rerun stage, electrostatic interactions are computed under 100 mM monovalent salt concentration, since open3SPN2 by default uses 100 mM monovalent salt concentration. 

Note for B_curved, the small energy difference in DNA bond, angle, and dihedral is because the different ways open3SPN2 and OpenSMOG3SPN2 build the template DNA structure if the sequence is W-C paired. The template built by OpenSMOG3SPN2 should be more accurate. Specifically, the template built by open3SPN2 is temp_temp.pdb, and the template built by OpenSMOG3SPN2 is cg_dna_template.pdb. For W-C paired sequence, open3SPN2 uses the full sequence (two ssDNA sequence combined) to build the template, and read the parameter from chain A of the built template. 

This test shows the code for DNA type A, B, and B_curved should all be correct. 

