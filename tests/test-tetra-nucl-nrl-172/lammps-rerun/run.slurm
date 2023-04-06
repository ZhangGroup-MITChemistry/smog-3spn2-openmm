#!/bin/bash 
#SBATCH --job-name=rerun
#SBATCH --output=rerun-%j.out 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=sched_mit_binz,newnodes,sched_mit_hill
#SBATCH --time=12:00:00
#SBATCH --mem=4G
#SBATCH --exclude=node[394]
#SBATCH --export=ALL
#SBATCH --constraint=centos7

module purge
module load gcc
module load engaging/openmpi/1.8.8

# print the node name that the job uses
echo "job uses node $SLURMD_NODENAME"

# For plumed 2.4.6 MATHEVAL package
export PLUMED_USE_LEPTON=yes

lammpsdir="/home/xclin/bin/awsem/bin"

mpirun -np 1 ${lammpsdir}/lmp_openmpi_plumed246_3spn_2016_sbm_rigidtemper_tanh -partition 1x1 -in rerun.in

echo "simulation finished"
echo "job done"

