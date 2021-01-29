#!/bin/bash
#SBATCH --time=00:45:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -A Research_Project-T109423
#SBATCH -p pq
#SBATCH --job-name=vcf_AF
#SBATCH --error=vcf_AF.err.txt
#SBATCH --output=vcf_AF.out.txt
#SBATCH --export=All
#SBATCH -D .


# this script is run on your local version HPC, make sure to change header/paths accordingly

module load VCFtools/0.1.16-foss-2018b-Perl-5.28.0
in_vcf=data/FIBR_STAR_pop_SNP.gatk.bi.miss.maf.final.filtered.depth4.recode.vcf.gz
out_dir=data/allele_freq
pop_files=data/pop_files

### getting relatedness and heterozygosity
vcftools --gzvcf $in_vcf \
--keep $pop_files/C_list.txt \
--freq \
--out $out_dir/IC_AF

vcftools --gzvcf $in_vcf \
--keep $pop_files/T_list.txt \
--freq \
--out $out_dir/IT_AF

vcftools --gzvcf $in_vcf \
--keep $pop_files/LL_list.txt \
--freq \
--out $out_dir/ILL_AF

vcftools --gzvcf $in_vcf \
--keep $pop_files/UL_list.txt \
--freq \
--out $out_dir/IUL_AF

vcftools --gzvcf $in_vcf \
--keep $pop_files/GH_list.txt \
--freq \
--out $out_dir/GHP_AF
