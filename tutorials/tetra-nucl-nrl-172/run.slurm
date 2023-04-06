#!/bin/bash
#SBATCH --job-name=chrom
#SBATCH --output=slurm-output/chrom-%j.out
#SBATCH --nodes=1
#SBATCH --gres=gpu:volta:1
#SBATCH --mem=10GB
#SBATCH --time=12:00:00

module load cuda/10.0

python tetra_nucl_nrl_172.py

echo "job done"

