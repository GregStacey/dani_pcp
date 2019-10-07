#!/bin/bash
#SBATCH -o prince_dani.out
#SBATCH -e prince_dani.err
#SBATCH --mail-user=richard.greg.stacey@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --time=023:55:00
#SBATCH --mem=16000M
#SBATCH --array=1-14

module load gcc/7.3.0
module load nixpkgs/16.09
module load netcdf/4.6.1
module load r/3.6.0

Rscript prince_dani.R $SLURM_ARRAY_TASK_ID
