#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1 # nodes=number of nodes required. ppn=number of processors per node
#SBATCH -A Research_Project-T109423
#SBATCH -p pq
#SBATCH --job-name=ROH
#SBATCH --error=ROH.err.txt
#SBATCH --output=ROH.out.txt
#SBATCH --export=All
#SBATCH -D .

# this script is run on your HPC, so change header/paths accordingly
# this script obtains the runs of homozygosity in each individual per population
# It requires phased vcfs, so make sure to run phase_vcfs.sh first

### set the variables
DIR=data/phased_vcfs
OUT=data/ROH

### --homozyg-snp is calculated as per the formula in Purfield et al 2017
~/lustre/bin/plink --vcf $DIR/GH_phased_autosomes.recode.vcf.gz --homozyg-snp 100 --homozyg-kb 500 --homozyg-density 50 --homozyg-gap 1000 --homozyg-window-snp 50 --homozyg-window-het 1 --homozyg-window-missing 5 --homozyg-window-threshold 0.05 --out $OUT/GHP_500K_1het --allow-extra-chr
~/lustre/bin/plink --vcf $DIR/LL_phased_autosomes.recode.vcf.gz --homozyg-snp 92 --homozyg-kb 500 --homozyg-density 50 --homozyg-gap 1000 --homozyg-window-snp 50 --homozyg-window-het 1 --homozyg-window-missing 5 --homozyg-window-threshold 0.05 --out $OUT/ILL_500K_1het --allow-extra-chr
~/lustre/bin/plink --vcf $DIR/UL_phased_autosomes.recode.vcf.gz --homozyg-snp 101 --homozyg-kb 500 --homozyg-density 50 --homozyg-gap 1000 --homozyg-window-snp 50 --homozyg-window-het 1 --homozyg-window-missing 5 --homozyg-window-threshold 0.05 --out $OUT/IUL_500K_1het --allow-extra-chr
~/lustre/bin/plink --vcf $DIR/C_phased_autosomes.recode.vcf.gz --homozyg-snp 131 --homozyg-kb 500 --homozyg-density 50 --homozyg-gap 1000 --homozyg-window-snp 50 --homozyg-window-het 1 --homozyg-window-missing 5 --homozyg-window-threshold 0.05 --out $OUT/C_500K_1het --allow-extra-chr
~/lustre/bin/plink --vcf $DIR/T_phased_autosomes.recode.vcf.gz --homozyg-snp 122 --homozyg-kb 500 --homozyg-density 50 --homozyg-gap 1000 --homozyg-window-snp 50 --homozyg-window-het 1 --homozyg-window-missing 5 --homozyg-window-threshold 0.05 --out $OUT/T_500K_1het --allow-extra-chr
