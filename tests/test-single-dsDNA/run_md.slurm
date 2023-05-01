#!/bin/bash
#SBATCH --job-name=md
#SBATCH --output=slurm-output/md-%j.out
#SBATCH --nodes=1
#SBATCH --gres=gpu:volta:1
#SBATCH --mem=20GB
#SBATCH --time=2-00:00:00

module load cuda/10.0

python run_md.py

echo "job done"

