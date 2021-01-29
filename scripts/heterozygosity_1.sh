#!/bin/bash
#SBATCH --time=00:45:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -A Research_Project-T109423
#SBATCH -p pq
#SBATCH --job-name=vcf_hwe_OI
#SBATCH --error=vcf_hwe_OI.err.txt
#SBATCH --output=vcf_hwe_OI.out.txt
#SBATCH --export=All
#SBATCH -D .

# this script is run on your HPC, so change header/paths accordingly

module load VCFtools/0.1.16-foss-2018b-Perl-5.28.0
in_vcf=data/FIBR_STAR_pop_SNP.gatk.bi.miss.maf.final.filtered.depth4.recode.vcf.gz
out_dir=data/heterozygosity
pop_files=data/pop_files

### getting relatedness and heterozygosity
vcftools --gzvcf $in_vcf \
--keep $pop_files/C_list.txt \
--hardy \
--out $out_dir/IC_het

vcftools --gzvcf $in_vcf \
--keep $pop_files/T_list.txt \
--hardy \
--out $out_dir/IT_het

vcftools --gzvcf $in_vcf \
--keep $pop_files/LL_list.txt \
--hardy \
--out $out_dir/ILL_het

vcftools --gzvcf $in_vcf \
--keep $pop_files/UL_list.txt \
--hardy \
--out $out_dir/IUL_het

vcftools --gzvcf $in_vcf \
--keep pop_files/GH_list.txt \
--hardy \
--out $out_dir/GHP_het
