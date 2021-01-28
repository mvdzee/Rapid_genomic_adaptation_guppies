# FIBR
To accompany FIBR paper

This repository contains various scripts that were used for population genomics analyses in  
**Rapid genomic adaptation in recently introduced populations of Trinidadian guppies (Poecilia reticulata)**

Raw reads were processed with a set of in-house scripts published in Fraser et al. (2020) - Improved reference genome uncovers novel sex-linked regions in the guppy (Poecilia reticulata)

All \*.sh scripts are written for SLURM based systems. Please change headers and batch format where necessary.

For selscan and Runs of Homozygosity, VCFs were phased using the following script:  
```phase_vcf.sh```  
For this you will need the final filtered bam files.  

Download the repo, it should have the following folders:
 1. scripts (this is where the scripts go)  

 2. data  
  - popgenome  
  - heterozygosity  
  - allele_freq  
  - ROH  
  - PCA  
  - selscan  
    - xpehh  
    - ihh12  
  - AF-vapeR  

 3. output  
  - popgenome  
  - heterozygosity  
  - allele_freq  
  - ROH  
  - PCA  
  - selscan  
    - xpehh  
    - ihh12 
  - AF-vapeR  
    
 4. figures  
  - popgenome  
  - heterozygosity  
  - allele_freq  
  - ROH  
  - PCA  
  - selscan  
    - xpehh  
    - ihh12 
  - AF-vapeR  

## Processing per analysis  
### Popgenome  
**To obtain raw data**  
Run these scripts on the HPC to get Fst, nuc.div. and Tajima's D:  
FST - ```popgenome.R```  
Tajima'S - ```popgenome.R```  
Nucleotide diversity (π) - ```popgenome.R```  
These scripts can be run by this shell:  ```popgenome.sh```  

**To process the raw data**  
On the HPC, in the folder with the outputs, run the following to concatenate the files per statistic:  
```
for stat in fst pi td; 
   do     
   awk 'FNR==1 && NR!=1 { while (/^chrom+/) getline; } 1 {print}' *${stat}*out > FIBR_${stat}.txt;
   done
```  
 Now move the output files to the main folder data/popgenome  
 
 To get FST means and medians run: ```fst_global.R```  
 To process Tajima's D run: ```TD_windows.R```  
 To process π run: ```Pi_windows.R```  

### VCFtools  
**To obtain raw data**  
Run these scripts on the HPC:  
Heterozygosity - ```heterozygosity_vcftools.sh```  
Allele frequency - ```allele_freq_vcftools.sh```  

**To process the data**  
*Heterozygosity*  
On the HPC, in the folder with the raw output, run the following:  
```
for POP in GHP GLP IC IT ILL IUL
   do
   cat ${POP}_het.hwe | tr ‘/’ ‘\t’ > ${POP}_het.txt
   done
```  
Now move the output files to the main folder data/heterozygosity 
and then run: ```Het_windows.R```  

*Plot genome-wide stats*  
To create the plot seen in fig 1D with values of expected heterozygosity, nucleotide diversity and Tajima's D, run  
``` plot_genomewide_stats.R ```  

*Allele frequency*  
On the HPC, in the folder with the raw output, run the following:  
```
for POP in GHP GLP IC IT ILL IUL
   do
   cat ${POP}_AF.frq | tr ':' '\t' | sed '1s/^.*$/CHROM\tPOS\tN_ALLELES\tN_CHR\tREF\tREF_FRQ\tALT\tALT_FRQ/' > ${POP}_AF.txt
   done
```  
Now move the \*.txt files to the main folder data/allele_freq  
and then run: ```AF_windows.R```  to calculate (absolute) delta allele frequencies.  

To calculate fixed/unchanged AFs between GHP and the LP populations run: ```fixed_AFs.R```  


### PCA
**To obtain the raw data**  
Run this scripts on the HPC:  
This script also prunes the VCF for linkage.  
```plink_PCA.sh```  

**To plot the PCA**  
Move the eigenvector and eigenvalue files from the HPC to the local data/PCA
Then run ```PCA_plot.R``` to create the plot in figure 1B


### Runs of homozygosity
This analysis requires phased VCFs!!  

**To obtain the raw data**  
Run these scripts on the HPC to get ROH per individual per population:  
Runs of homozygosity - ```plink_ROH.sh```  

**To process the raw data**  
Move the output files from the HPC output folder to the local FIBR_ms/data/ROH
To process the files and create the plot seen in figure 1C run: ```ROH_plots.R``` 


### Selscan genome scans
This analysis requires phased VCFs!!  

**To obtain the raw data**  
Run this scripts on the HPC to get the required plink \*.map files in the correct format as well as XP-EHH and iHH12:  
Selscan - ```selscan.sh```  

**To process the raw data**  
*XP-EHH*  
On the HPC, in the folder with the XP-EHH outputs, run the following to concatenate the files per population:  
```
for POP in C T LL UL; 
   do     
   awk 'FNR==1 && NR!=1 { while (/^id+/) getline; } 1 {print}' GH_${POP}*out.norm > GH_${POP}_XPEHH.txt;
   done
 ```  
 
Move the concatenated file from the HPC to the local data/selscan/xpehh folder.
In R, run ```XPEHH_outliers.R```  
This script will output the outliers per introduced population.  
 
*iHH12*  
On the HPC, in the folder with the iHH12 outputs, run the following to concatenate the files per population:  
```
for POP in GH C T LL UL; 
   do     
   awk 'FNR==1 && NR!=1 { while (/^id+/) getline; } 1 {print}' ${POP}*out.norm > ${POP}_iHH12.txt;
   done
 ```  
 
Move the concatenated file from the HPC to the local data/selscan/ihh12 folder.
In R, run ```ihh12_outliers.R```  
This script will output the outliers per introduced population.

To plot the selscan results:
Run the following script: ```selscan_plots.R```  
This will create the plots seen in figure 2.  


### AF-vapeR  
There will a more detailed description of this software available on \*jims github\*.  
**To obtain the raw data and process**  
For this analysis you will need the scripts:  
```AF_vapeR_chrom.R``` and ```AF_vapeR_functions.R``` in the same folder   
You will also need to have BCFtools installed. For a version without the use of BCFtools please see \*jims github\*.  

In RStudio, run through the ```AF_vapeR_chrom.R``` to get the results and figures.  
