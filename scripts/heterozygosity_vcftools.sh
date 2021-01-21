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
in_vcf=/gpfs/ts0/home/mv323/lustre/start_up_data/FIBR/STAR/data/FIBR_gvcfs/FIBR_STAR_pop_SNP.gatk.bi.miss.maf.final.filtered.depth4.recode.vcf
out_dir=/gpfs/ts0/home/mv323/lustre/start_up_data/people/mijke/STAR_analyses/FIBR/heterozygosity

### getting relatedness and heterozygosity
vcftools --vcf $in_vcf \
--keep /gpfs/ts0/home/mv323/lustre/scripts/vcftools/c_list.txt \
--hardy \
--out $out_dir/IC_het

vcftools --vcf $in_vcf \
--keep /gpfs/ts0/home/mv323/lustre/scripts/vcftools/t_list.txt \
--hardy \
--out $out_dir/IT_het

vcftools --vcf $in_vcf \
--keep /gpfs/ts0/home/mv323/lustre/scripts/vcftools/ll_list.txt \
--hardy \
--out $out_dir/ILL_het

vcftools --vcf $in_vcf \
--keep /gpfs/ts0/home/mv323/lustre/scripts/vcftools/ul_list.txt \
--hardy \
--out $out_dir/IUL_het

vcftools --vcf $in_vcf \
--keep /gpfs/ts0/home/mv323/lustre/scripts/vcftools/gl_list.txt \
--hardy \
--out $out_dir/GLP_het

vcftools --vcf $in_vcf \
--keep /gpfs/ts0/home/mv323/lustre/scripts/vcftools/gh_list.txt \
--hardy \
--out $out_dir/GHP_het