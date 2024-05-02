#!/bin/bash
#SBATCH --time=144:00:00
#SBATCH --partition=week-long
#SBATCH --nodes=1
#SBATCH --mem=1G
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpuq
echo $HOSTNAME
echo $SLURM_SUBMIT_DIR
echo $SLURM_SUBMIT_HOST
echo $SLURM_JOB_ID
autosolvate mdrun -f Fe_plus2_solvated
#module load orca/5.0.2
##python FFmetalcomplex.py -n Fe_plus2 -u 1 -c 2 -x orca -p 32 -G /opt/orca/5.0.2/orca -e acetonitrile
