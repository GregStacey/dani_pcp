#!/bin/bash
#SBATCH -o prince_dani.out
#SBATCH -e prince_dani.err
#SBATCH --mail-user=richard.greg.stacey@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --time=02:55:00
#SBATCH --mem=16000M

module load gcc/7.3.0
module load nixpkgs/16.09
module load netcdf/4.6.1
module load r/3.6.0

PROJECT_DIR=~/projects/def-ljfoster/rstacey/dani_pcp/
cd ${PROJECT_DIR}/R

Rscript prince_dani.R
