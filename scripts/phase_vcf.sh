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
# install beagle and shapeit following their instructions
# Set up the environment
module load Java/1.8.0_74
module load VCFtools/0.1.16-foss-2018b-Perl-5.28.0 BCFtools/1.9-foss-2018b
module load HTSlib/1.8-foss-2018a

MASTER=data/phased_vcfs
POPMAP_DIR=data/pop_files
DATASET=FIBR
VCF=data/FIBR_finalVCF.vcf.gz

# a list of the chromosomes/scaffolds to process
samples=$MASTER/chroms.txt
chrs=(`cat $samples`)
chr_array=${chrs[(($SLURM_ARRAY_TASK_ID))]}

# array with your population names
pop_array=(GH T LL UL C)
for POP in "${pop_array[@]}"
do

# change the path to a folder that contains a sub-folder per population with the bam files for that population
BAM_DIR=path/to/bamfiles   
# Extract
vcftools --gzvcf $VCF \
--keep $POPMAP_DIR/${POP}_list.txt \
--chr $chr_array \
--recode \
--out $MASTER/${POP}_${chr_array}_tmp

#-----------------------------
# Phase with beagle
java -Xss20m -jar ~/lustre/bin/beagle.28Sep18.793.jar \     # change the path to where beagle is
gt=$MASTER/${POP}_${chr_array}_tmp.recode.vcf \
out=$MASYER/${POP}_${chr_array}_beagle_phased.gt \
nthreads=12

#Tabix index
tabix $MASTER/${POP}_${chr_array}_beagle_phased.gt.vcf.gz

#-----------------------------
#Phase with shapeit
#First we need to list the BAMs
rm -f $MASTER/${POP}_${CHR}_bamlist
for ind in $(cat $POPMAP_DIR/${POP}_list.txt)
do
ind2=$(echo -e ${ind}_)
bam_in=$(ls $BAM_DIR | grep "${ind2}" | grep ".bam$")
echo -e "${ind}\t${bam_in}\t$chr_array" >> $MASTER/${POP}_${chr_array}_bamlist
done

BAM_LIST=$MASTER/${POP}_${chr_array}_bamlist

cd $BAM_DIR

/gpfs/ts0/home/mv323/lustre/software/extractPIRs.v1.r68.x86_64/extractPIRs --bam $BAM_LIST \      # change the path to the correct folder for extractPIRs
           --vcf $MASTER/${POP}_${chr_array}_beagle_phased.gt.vcf.gz \
           --out $MASTER/${POP}_${chr_array}_beagle_PIRlist

## Now phase
/gpfs/ts0/home/mv323/lustre/software/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit -assemble \     # change the path to the correct folder for shapeit
--force \
        --input-vcf $MASTER/${POP}_${chr_array}_beagle_phased.gt.vcf.gz \
        --input-pir $MASTER/${POP}_${chr_array}_beagle_PIRlist \
        -O $MASTER/${POP}_${chr_array}_shapeit_beagle_haps

## And convert to VCF
/gpfs/ts0/home/mv323/lustre/software/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit -convert \     # change the path to the correct folder for shapeit
        --input-haps $MASTER/${POP}_${chr_array}_shapeit_beagle_haps \
        --output-vcf $MASTER/${POP}_${chr_array}_shapeit_beagle.vcf

bgzip $MASTER/${POP}_${chr_array}_shapeit_beagle.vcf > $MASTER/${POP}_${chr_array}_shapeit_beagle.vcf.gz
tabix -p vcf $MASTER/${POP}_${chr_array}_shapeit_beagle.vcf.gz

done
