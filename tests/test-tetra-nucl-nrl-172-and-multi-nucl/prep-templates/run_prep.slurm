#!/bin/bash
#SBATCH --job-name=prep
#SBATCH --output=slurm-output/prep-%j.out
#SBATCH --nodes=1
#SBATCH --mem=10GB
#SBATCH --time=12:00:00

python prep_x3dna_template.py

echo "job done"

