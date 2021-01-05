#!/bin/bash
#SBATCH --time=01:0:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1 # nodes=number of nodes required. ppn=number of processors per node
#SBATCH -A Research_Project-T109423
#SBATCH -p pq
#SBATCH --job-name=pca_FIBR
#SBATCH --error=pca_FIBR.err.txt
#SBATCH --output=pca_FIBR.out.txt
#SBATCH --export=All
#SBATCH -D .

# this script is run on your HPC, so change header/paths accordingly

# set variables & paths
vcfs=~/lustre/start_up_data/FIBR/STAR/data/FIBR_gvcfs
out=/gpfs/ts0/home/mv323/lustre/start_up_data/people/mijke/STAR_analyses/FIBR/plink/PCA

#### Get the lists of sites that are in LE (prune.in)
/gpfs/ts0/home/mv323/lustre/bin/plink --vcf $vcfs/FIBR_STAR_pop_SNP.gatk.bi.miss.maf.final.filtered.depth4.recode.vcf.gz \
--allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 5 0.2 \
--out $out/FIBR_PCA

### Extract the LE sits from the oringal vcf
/gpfs/ts0/home/mv323/lustre/bin/plink --vcf $vcfs/FIBR_STAR_pop_SNP.gatk.bi.miss.maf.final.filtered.depth4.recode.vcf.gz \
--extract $out/FIBR_PCA.prune.in \
--allow-extra-chr \
--set-missing-var-ids @:# \
--recode vcf \
--out $vcfs/pruned/FIBR_PCA_pruned

### calculate eigenvecs & vals for the pruned vcf
/gpfs/ts0/home/mv323/lustre/bin/plink --vcf $vcfs/pruned/FIBR_PCA_pruned.vcf \
--allow-extra-chr \
--pca header tabs \
--out $out/FIBR_PCA
