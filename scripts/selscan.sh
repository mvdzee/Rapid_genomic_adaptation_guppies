#!/bin/bash
#SBATCH --time=40:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4 # nodes=number of nodes required. ppn=number of processors per node
#SBATCH -A Research_Project-T109423
#SBATCH -p pq
#SBATCH --job-name=selscan
#SBATCH --error=selscan.err.txt
#SBATCH --output=selscan.out.txt
#SBATCH --export=All
#SBATCH -D .

# Set up the environment
VCFDIR=/gpfs/ts0/home/mv323/lustre/start_up_data/people/mijke/STAR_analyses/FIBR/phased_vcfs/outputs
OUT=/gpfs/ts0/home/mv323/lustre/start_up_data/people/mijke/STAR_analyses/FIBR/selscan
samples=$OUT/chroms.txt
pop_array=(C T LL UL GH GL)

for POP in "${pop_array[@]}"
do
arr=()
while IFS= read -r line;
do
    arr+=("$line")
done <  $samples

for i in ${arr[@]}
do
~/lustre/bin/plink --vcf $VCFDIR/chroms/${POP}_${i}_shapeit_beagle.vcf.gz --recode --out $VCFDIR/chroms/${POP}_${i}_shapeit_beagle.vcf.gz
cat $VCFDIR/chroms/${POP}_${i}_shapeit_beagle.vcf.gz.map | awk -v var1="$i" -v var2="$POP" '{print var1"\t"var1":"$4"\t"$4"\t"$4}' > test && mv test $VCFDIR/chroms/${POP}_${i}_shapeit_beagle.vcf.gz.map
done

# calculate iHH12
for i in ${arr[@]}
do
~/lustre/bin/selscan --ihh12 --vcf $VCFDIR/chroms/${POP}_${i}_shapeit_beagle.vcf.gz --map $VCFDIR/chroms/${POP}_${i}_shapeit_beagle.vcf.gz.map --maf 0.01 --trunc-ok --out $OUT/ihh12/${POP}_${i}_ihh12 --threads 4
done

# calculate XP-EHH
for i in ${arr[@]}
do
~/lustre/bin/selscan --xpehh --vcf $VCFDIR/chroms/C_${i}_shapeit_beagle.vcf.gz --vcf-ref $VCFDIR/chroms/GH_${i}_shapeit_beagle.vcf.gz --map $VCFDIR/chroms/GH_${i}_shapeit_beagle.vcf.gz.map --trunc-ok --maf 0.01 --out $OUT/xpehh/GH_C_${i} --threads 4
~/lustre/bin/selscan --xpehh --vcf $VCFDIR/chroms/T_${i}_shapeit_beagle.vcf.gz --vcf-ref $VCFDIR/chroms/GH_${i}_shapeit_beagle.vcf.gz --map $VCFDIR/chroms/GH_${i}_shapeit_beagle.vcf.gz.map --trunc-ok --maf 0.01 --out $OUT/xpehh/GH_T_${i} --threads 4
~/lustre/bin/selscan --xpehh --vcf $VCFDIR/chroms/LL_${i}_shapeit_beagle.vcf.gz --vcf-ref $VCFDIR/chroms/GH_${i}_shapeit_beagle.vcf.gz --map $VCFDIR/chroms/GH_${i}_shapeit_beagle.vcf.gz.map --trunc-ok --maf 0.01 --out $OUT/xpehh/GH_LL_${i} --threads 4
~/lustre/bin/selscan --xpehh --vcf $VCFDIR/chroms/UL_${i}_shapeit_beagle.vcf.gz --vcf-ref $VCFDIR/chroms/GH_${i}_shapeit_beagle.vcf.gz --map $VCFDIR/chroms/GH_${i}_shapeit_beagle.vcf.gz.map --trunc-ok --maf 0.01 --out $OUT/xpehh/GH_UL_${i} --threads 4
~/lustre/bin/selscan --xpehh --vcf $VCFDIR/chroms/GL_${i}_shapeit_beagle.vcf.gz --vcf-ref $VCFDIR/chroms/GH_${i}_shapeit_beagle.vcf.gz --map $VCFDIR/chroms/GH_${i}_shapeit_beagle.vcf.gz.map --trunc-ok --maf 0.01 --out $OUT/xpehh/GH_GL_${i} --threads 4
done
done

# to normalise, for each pop do:
for POP in C T LL UL GL
do
~/lustre/software/selscan-master/bin/linux/norm --xpehh --crit-percent 0.05 --winsize 75000 --files $VCFDIR/xpehh/GH_${POP}*out
done

for POP in C T LL UL GH GL
do
~/lustre/software/selscan-master/bin/linux/norm --ihh12 --crit-percent 0.05 --winsize 75000 --files $VCFDIR/ihh12/${POP}*out
done