#!/bin/bash
#SBATCH --job-name=nucl
#SBATCH --output=slurm-output/nucl-%j.out
#SBATCH --nodes=1
#SBATCH --gres=gpu:volta:1
#SBATCH --mem=10GB
#SBATCH --time=12:00:00

module load cuda/10.0

python run.py

echo "job done"

