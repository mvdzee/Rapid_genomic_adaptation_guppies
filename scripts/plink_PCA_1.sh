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


# this script prunes the vcf for linkage and calculates the eigenvalues/vectors
VCF=data/FIBR_finalVCF.vcf.gz
OUT=data/PCA

#### Get the lists of sites that are in LE (prune.in)
/gpfs/ts0/home/mv323/lustre/bin/plink --vcf $VCF \
--allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 5 0.2 \
--out $OUT/FIBR_PCA

### Extract the LE sits from the oringal vcf
/gpfs/ts0/home/mv323/lustre/bin/plink --vcf $VCF \
--extract $OUT/FIBR_PCA.prune.in \
--allow-extra-chr \
--set-missing-var-ids @:# \
--recode vcf \
--out $OUT/FIBR_PCA_pruned

# ### calculate eigenvecs & vals for the pruned vcf
/gpfs/ts0/home/mv323/lustre/bin/plink --vcf $OUT/FIBR_PCA_pruned.vcf \
--allow-extra-chr \
--pca header tabs \
--out $OUT/FIBR_PCA
