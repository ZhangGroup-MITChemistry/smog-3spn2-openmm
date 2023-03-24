#!/usr/bin/env python

import os
import sys
import shutil
from subprocess import *
import fileinput
from numpy import *

################################################

def my_lt_range(start, end, step):
    while start < end:
        yield start
        start += step

def my_le_range(start, end, step):
    while start <= end:
        yield start
        start += step

#############################################

globalDir = os.getcwd() # get current directory

# Add line to the top of a file;
def line_prepender(filename, line):
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0) # locate to the beginning of f
        f.write(line.rstrip('\r\n') + '\n' + content)

def create_folder(runId): # create folder according to runId and put plumed script into it 
    # prepare folders
    subdir = "sim" 
    simFolder = "%s/%s/" % (globalDir, subdir)
    rundir = simFolder + "run%03d/" % runId 

    if os.path.exists(rundir):
        shutil.rmtree(rundir) 
        os.makedirs(rundir)
    else:
        os.makedirs(rundir)

    # create plumed script
    fh = open(rundir + 'plumed.txt', 'w')
    cv_tmp = fileinput.input('plumed.txt')

    for line in cv_tmp:
        fh.write(line)
    fh.close()

    # If it is a restarting run, append "RESTART" in the first line, required by the PLUMED format
    if (runId != 0):
        line_prepender(rundir + 'plumed.txt', 'RESTART')

def submit(staId, endId):
    subdir = "sim"
    simFolder = "%s/%s/" % (globalDir, subdir)

    # create pbs file
    for i in range(staId, endId + 1):
        pbsFile = simFolder + "myjob_eofe-%d.slurm" % i
        pbs_tmp = fileinput.input('myjob_eofe.slurm')
        pf = open(pbsFile, 'w')

        # add head to myjob_eofe.slurm
        head = '#!/bin/bash \n#SBATCH --job-name=tetranucl \n#SBATCH --output=slurm-%j-curr-' + str(i) + '.out \n'
        pf.write(head)

        # set curr and max_id
        for line in pbs_tmp:
            if '#INPUT_CURR_AND_MAX_ID' in line:
                pf.write('curr=%d \n' % i)
                pf.write('max_id=%d \n' % endId)
            else:
                pf.write(line)
        pf.close()

    # submit the job 
    cmd = "cd %s; sbatch myjob_eofe-%d.slurm; cd %s" % (simFolder, staId, globalDir) 
    q = Popen(cmd, shell=True, stdout=PIPE) 
    q.communicate()

if __name__ == '__main__':
    nsta = int(sys.argv[1])
    nend = int(sys.argv[2])
    for runId in my_le_range(nsta, nend, 1):
        create_folder(runId)

    submit(nsta, nend)

