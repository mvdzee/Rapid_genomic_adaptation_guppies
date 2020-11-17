# FIBR
To accompany FIBR paper

This repository contains various scripts that were used for population genomics analyses in  
**Title of paper**

Raw reads were processed with a set of in-house scripts published in Fraser et al. (2020) - Improved reference genome uncovers novel sex-linked regions in the guppy (Poecilia reticulata)

Download the scripts into a folder (FIBR_ms) with the following sub-folders: 
 1. scripts (this is where the scripts go)  

 2. data  
  - popgenome  
  - heterozygosity  
  - allele_freq  
  - selscan  
    - xpehh  
    - ihh12  

 3. output  
  - popgenome  
  - heterozygosity  
  - allele_freq  
  - selscan  
    - xpehh  
    - ihh12 
    
 4. figures  

### To obtain raw data for all analyses  
These scripts are run on the HPC where each analysis has its own folder (so one for popgenome, one for heterozygosity, one for allele frequency and one for selscan)
FST - ```popgenome.R```  
Tajima'S - ```popgenome.R```  
Nucleotide diversity (π) - ```popgenome.R```  
Heterozygosity - ```heterozygosity_vcftools.sh```  
Allele frequency - ```allele_freq_vcftools.sh```  
Selscan XP-EHH -  
Selscan iHH12 - 

### To process raw data  

#### Popgenome results   
On the HPC, in the folder with the outputs, run the following to concatenate the files per statistic:  
```for stat in fst pi td; 
   do     
   for file in *${stat}*out;
      do
      awk 'FNR==1 && NR!=1 { while (/^chrom+/) getline; } 1 {print}' *${stat}*out > FIBR_${stat}.txt;
      done;
   done
   ```  
 Now move the output files to the main folder FIBR_ms/data/popgenome  
 
 To get FST means and medians run: ```fst_global.R```  
 To process Tajima's D run: ```TD_windows.R```  
 To process π run: ```Pi_windows.R```  
 
#### VCFtools results  
*Heterozygosity*  
On the HPC, in the folder with the raw output, run the following:  
```for POP in GHP GLP IC IT ILL IUL
   do
   cat ${POP}_het.hwe | tr ‘/’ ‘\t’ > ${POP}_het.txt
   done
```  
Now move the output files to the main folder FIBR_ms/data/heterozygosity 
and then run: ```Het_windows.R```  

*Allele frequency*  
On the HPC, in the folder with the raw output, run the following:  
```for POP in GHP GLP IC IT ILL IUL
   do
   cat ${POP}_AF.frq | tr ':' '\t' | sed '1s/^.*$/CHROM\tPOS\tN_ALLELES\tN_CHR\tREF\tREF_FRQ\tALT\tALT_FRQ/' > ${POP}_AF.txt
   done
```  
Now move the output files to the main folder FIBR_ms/data/allele_freq  
and then run: ```AF_windows.R```  to calculate (absolute) delta allele frequencies.  

To calculate fixed/unchanged AFs between GHP and the LP populations run: ```fixed_AFs.R```  

#### Plot neutral stats 
