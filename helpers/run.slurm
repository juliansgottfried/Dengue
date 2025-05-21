#!/bin/bash
#SBATCH --partition=cs
#SBATCH --array=1-50
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem-per-cpu=10GB
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=juliansgottfried@gmail.com

module unload r
module unload gcc

module load r/gcc/4.3.1
Rscript helpers/fit.R \
    $SLURM_ARRAY_TASK_MAX \
    $SLURM_ARRAY_TASK_ID \
    $SLURM_NTASKS_PER_NODE \
    $SLURM_JOB_NAME
