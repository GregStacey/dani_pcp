#!/bin/bash
#SBATCH -o experiment-predict-reproducibility-ms.out
#SBATCH -e experiment-predict-reproducibility-ms.err
#SBATCH --mail-user=richard.greg.stacey@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --time=023:55:00
#SBATCH --mem=16000M
#SBATCH --array=1-14

module load r/3.4.3
module load gcc/5.4.0

Rscript prince_dani.R $SLURM_ARRAY_TASK_ID