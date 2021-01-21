#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH -A Research_Project-T109423
#SBATCH -p pq
#SBATCH --job-name=phasing
#SBATCH --error=phasing.err.txt
#SBATCH --output=phasing.out.txt
#SBATCH --export=All
#SBATCH -D .
#SBATCH --array=0-266%20

# this script phases your VCF files, first using beagle and then a second round with shapeit
# Set up the environment
module load Java/1.8.0_74
module load VCFtools/0.1.16-foss-2018b-Perl-5.28.0 BCFtools/1.9-foss-2018b
module load HTSlib/1.8-foss-2018a

MASTER=/gpfs/ts0/home/mv323/lustre/start_up_data/people/mijke/STAR_analyses/FIBR/phased_vcfs
OUT=$MASTER/tmp
POPMAP_DIR=/gpfs/ts0/home/mv323/lustre/popmap/FIBR
DATASET=FIBR
VCF=/gpfs/ts0/home/mv323/lustre/start_up_data/FIBR/STAR/data/FIBR_gvcfs/FIBR_STAR_pop_SNP.gatk.bi.miss.maf.final.filtered.depth4.recode.vcf.gz



# a list of the chromosomes/scaffolds to process
samples=$MASTER/scafs.txt
chrs=(`cat $samples`)
chr_array=${chrs[(($SLURM_ARRAY_TASK_ID))]}

# array with your population names
pop_array=(GH GL T LL UL C)
for POP in "${pop_array[@]}"
do

BAM_DIR=/gpfs/ts0/projects/Research_Project-T110748/STAR_bams
# Extract
vcftools --gzvcf $VCF \
--keep $POPMAP_DIR/${POP}.popmap \
--chr $chr_array \
--recode \
--out $OUT/${POP}_${chr_array}_tmp

#-----------------------------
# Phase with beagle
java -Xss20m -jar ~/lustre/bin/beagle.28Sep18.793.jar \
gt=$OUT/${POP}_${chr_array}_tmp.recode.vcf \
out=$OUT/tmp2/${POP}_${chr_array}_beagle_phased.gt \
nthreads=12

#Tabix index
tabix $OUT/tmp2/${POP}_${chr_array}_beagle_phased.gt.vcf.gz

#-----------------------------
#Phase with shapeit
#First we need to list the BAMs
rm -f $MASTER/outputs/${POP}_${CHR}_bamlist
for ind in $(cat $POPMAP_DIR/${POP}.popmap)
do
ind2=$(echo -e ${ind}_)
bam_in=$(ls $BAM_DIR | grep "${ind2}" | grep ".bam$")
echo -e "${ind}\t${bam_in}\t$chr_array" >> $MASTER/outputs/${POP}_${chr_array}_bamlist
done

BAM_LIST=$MASTER/outputs/${POP}_${chr_array}_bamlist

cd $BAM_DIR

/gpfs/ts0/home/mv323/lustre/software/extractPIRs.v1.r68.x86_64/extractPIRs --bam $BAM_LIST \
           --vcf $MASTER/tmp/tmp2/${POP}_${chr_array}_beagle_phased.gt.vcf.gz \
           --out $OUT/${POP}_${chr_array}_beagle_PIRlist

## Now phase
/gpfs/ts0/home/mv323/lustre/software/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit -assemble \
--force \
        --input-vcf $MASTER/tmp/tmp2/${POP}_${chr_array}_beagle_phased.gt.vcf.gz \
        --input-pir $OUT/${POP}_${chr_array}_beagle_PIRlist \
        -O $MASTER/outputs/intermediate/${POP}_${chr_array}_shapeit_beagle_haps

## And convert to VCF
/gpfs/ts0/home/mv323/lustre/software/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit -convert \
        --input-haps $MASTER/outputs/intermediate/${POP}_${chr_array}_shapeit_beagle_haps \
        --output-vcf $MASTER/outputs/${POP}_${chr_array}_shapeit_beagle.vcf

bgzip $MASTER/outputs/${POP}_${chr_array}_shapeit_beagle.vcf > $OUT/${POP}_${chr_array}_shapeit_beagle.vcf.gz
tabix -p vcf $OUT/${POP}_${chr_array}_shapeit_beagle.vcf.gz
# #
# #cd $MASTER
# #
# ## Tidy
# mv $OUT/${POP}_${chr_array}_*vcf* $MASTER/outputs/
#rm -f $OUT/${POP}_${chr_array}_*
done