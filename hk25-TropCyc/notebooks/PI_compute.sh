#!/bin/bash
#SBATCH -J PI
#SBATCH --time=24:00:00
#SBATCH --ntasks=16
#SBATCH --no-requeue
#SBATCH --mem=64G
#SBATCH --export=NONE
#SBATCH -N 1
#SBATCH --sockets-per-node=1
#SBATCH --account=bb1153

#unset SLURM_EXPORT_ENV

#module load python3/2023.01-gcc-11.2.0
module load openmpi/4.1.6-oneapi-2024.2.1
source /home/b/b383007/.bashrc
source activate hk25
#conda activate /home/b/b383007/.conda/envs/hk25
python3 --version
which python3

mpirun -np 16 --bind-to core:overload-allowed python3 PI_compute.py 

