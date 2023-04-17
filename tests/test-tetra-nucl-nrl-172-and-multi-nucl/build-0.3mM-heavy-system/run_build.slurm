#!/bin/bash
#SBATCH --job-name=build
#SBATCH --output=slurm-output/build-%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=40GB
#SBATCH --time=12:00:00

python build_nonrigid_heavy_system.py

echo "job done"

