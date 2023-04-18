#!/bin/bash
#SBATCH --job-name=md
#SBATCH --output=slurm-output/md-%j.out
#SBATCH --nodes=1
#SBATCH --gres=gpu:volta:1
#SBATCH --mem=10GB
#SBATCH --time=2-00:00:00

module load cuda/10.0

python step1_relax.py --output_dcd relax.dcd --output_interval 10000 --n_steps 100000

python step2_prod.py --input_dcd relax.dcd --output_dcd prod.dcd --output_interval 10000 --n_steps 1000000

echo "job done"

