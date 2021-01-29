#!/bin/bash
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -A Research_Project-T109423
#SBATCH -p pq
#SBATCH --job-name=vcf_hwe_OI
#SBATCH --error=popgenome.err.txt
#SBATCH --output=popgenome.out.txt
#SBATCH --export=All
#SBATCH -D .

# this script is run on your HPC, so change header/paths accordingly
# load R module
module load R/3.4.3-foss-2016b-X11-20160819

Rscript popgenome.R
