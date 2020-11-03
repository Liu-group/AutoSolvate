#!/bin/bash
#SBATCH --job-name="gauss"
#SBATCH --output="gauss.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --export=ALL
#SBATCH -t 00:30:00

source ~/.bashrc
export GAUSS_SCRDIR=/scratch/$USER/$SLURM_JOBID
conda activate autoneb
module load amber/20
module unload python/2.7.15
module load gnu
module load gaussian/16.C.01
python set_up_solvated_AmberMD.py 
