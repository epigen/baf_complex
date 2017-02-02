#!/bin/bash
#SBATCH --partition=shortq
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32000

date
module unload R
module load R/3.2.3

/cm/shared/apps/R/3.2.3/bin/Rscript ~/jobs/copywriteR.baf.R $1 $2

date
