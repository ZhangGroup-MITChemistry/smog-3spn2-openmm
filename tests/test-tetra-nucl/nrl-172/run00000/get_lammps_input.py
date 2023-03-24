#!/usr/bin/env python

import sys
from numpy import *
from subprocess import *
import fileinput

def genLammpsInput(restartFolder, outFolder, runId):
    lammps_template = '../../template.lammps' 
    
    # We need to decide whether it is the first run or a restarting run;
    if (runId == 0):
        fh = open(outFolder + "/proteinDna_sim.in", "w") 
        for line in open(lammps_template).readlines():
            items = line.split()
            fh.write(line)
            if len(items) >= 2 and items[0] == "read_data":
                fh.write("read_dump           ../../DUMP_initial.dcd 0 x y z box no format molfile dcd /home/xclin/lib/vmd.old/plugins/LINUXAMD64/molfile\n")
        fh.close()

    else: # for the restarting run, read the final snapshot of the former run
        fh = open(outFolder + "/proteinDna_sim.in", "w")
        for line in open(lammps_template).readlines():
            items = line.split()

            if len(items) >= 2 and items[0] == "read_data" and runId >= 1:
                fh.write(line)
                cmd = '/home/xclin/bin/anaconda2/bin/mdconvert -o DUMP_Extract.dcd -f -i %d ../run%03d/DUMP_FILE.dcd' % (-1, runId - 1) 
                q = Popen(cmd, shell=True, stdout=PIPE)
                q.communicate()
                fh.write("read_dump           DUMP_Extract.dcd 0 x y z format molfile dcd /home/xclin/lib/vmd.old/plugins/LINUXAMD64/molfile\n")
                
            elif runId >= 1 and len(items) >= 2 and items[0] == "minimize":
                pass
            else:
                fh.write(line)
    
if __name__ == "__main__":
    genLammpsInput(sys.argv[1], sys.argv[2], int(sys.argv[3]))

