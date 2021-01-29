### Rscript for popGenome Dxy ###
library("PopGenome", lib.loc="/gpfs/ts0/home/mv323/lustre/bin")
library(data.table)
library(dplyr)

#Import fai for genome data
chr_length<-read.table("data/popgenome/FIBR_chroms.fai",header=F)

####### Loop over chromosomes ###########
chr_function<-function(x){

#chr name as defined from fasta.fai
tid=as.character(chr_length[x,]$V1)

#chr length
topos = chr_length[x,]$V2

#Read in VCF
vcf_file <- readVCF('data/FIBR_finalVCF.vcf.gz', tid=tid, frompos = 1, topos = topos, numcols = 1000000, include.unknown=TRUE)

### assign pops ###
c_pop <- c('C1','C4','C6','C8','C9','C12','C13','C14','C15','C16','C17','C20','C21','C22','C23')
gh_pop <- c('GH1','GH2','GH3','GH4','GH5','GH6','GH7','GH8','GH9','GH10','GH11','GH12','GH13','GH14','GH15','GH17','GH18','GH19','GH20')
ll_pop <- c('LL1','LL2','LL3','LL5','LL6','LL7','LL9','LL10','LL11','LL13','LL14')
t_pop <- c('T2','T3','T4','T5','T6','T7','T8','T10','T11','T13','T14','T15','T16','T19')
ul_pop <- c('UL1','UL2','UL3','UL6','UL7','UL8','UL9','UL10','UL11','UL13','UL15','UL17','UL18','UL19','UL20')
vcf_file <- set.populations(vcf_file, list(c_pop,gh_pop,ll_pop,t_pop,ul_pop), diploid=TRUE)

### create windows ###
slide_vcf_file_75kb <- sliding.window.transform(vcf_file, 75000, 75000, type=2)

### do pop stats
slide_vcf_file_75kb <- F_ST.stats(slide_vcf_file_75kb, mode="nucleotide")
slide_vcf_file_75kb <- neutrality.stats(slide_vcf_file_75kb, FAST=TRUE)
slide_vcf_file_75kb <- diversity.stats.between(slide_vcf_file_75kb, nucleotide.mode=TRUE)
slide_vcf_file_75kb <- detail.stats(slide_vcf_file_75kb,site.spectrum=TRUE)

### create output for 75kb files ###
calc_FST<-function(fst_list){

        tmp<-data.frame(IC_GHP = fst_list[,1],
                IC_ILL = fst_list[,2],
                IC_IT = fst_list[,3],
                IC_IUL = fst_list[,4],
                GHP_ILL = fst_list[,5],
                GHP_IT = fst_list[,6],
                GHP_IUL = fst_list[,7],
                ILL_IT = fst_list[,8],
                ILL_IUL = fst_list[,9],
                IT_IUL = fst_list[,10])
        if(length(tmp) == 0){
        return('NA')
        } else {
        tmp['chrom']<-rep(tid,nrow(tmp))
        tmp['window']<-1:nrow(tmp)
        tmp['window_start']<- ((tmp$window-1)*75000)+1
        tmp['window_end']<- tmp$window*75000
        tmp<-tmp[,c(11:14,1:10)]
        }
        }
FST_mat<-data.frame(calc_FST(t(slide_vcf_file_75kb@nuc.F_ST.pairwise)))
outfile_fst <-paste0("data/popgenome/",tid,".fst.popgenome.out")
write.table(FST_mat, outfile_fst, sep = "\t", row.names=F, quote = F)

calc_PI<-function(pi_list){

        tmp<-data.frame(IC = pi_list[,1]/75000,
                GHP = pi_list[,2]/75000,
                ILL = pi_list[,3]/75000,
                IT = pi_list[,4]/75000,
                IUL = pi_list[,5]/75000
                )
        if(length(tmp) == 0){
        return('NA')
        } else {
        tmp['chrom']<-rep(tid,nrow(tmp))
        tmp['window']<-1:nrow(tmp)
        tmp['window_start']<- ((tmp$window-1)*75000)+1
        tmp['window_end']<- tmp$window*75000
        tmp<-tmp[,c(6:9,1:5)]
        }
        }
PI_mat<-data.frame(calc_PI(slide_vcf_file_75kb@nuc.diversity.within))
outfile_pi <-paste0("data/popgenome/",tid,".pi.popgenome.out")
write.table(PI_mat, outfile_pi, sep = "\t", row.names=F, quote = F)

calc_TD<-function(td_list){

        tmp<-data.frame(IC = td_list[,1],
                GHP = td_list[,2],
                ILL = td_list[,3],
                IT = td_list[,4],
                IUL = td_list[,5])
        if(length(tmp) == 0){
        return('NA')
        } else {
        tmp['chrom']<-rep(tid,nrow(tmp))
        tmp['window']<-1:nrow(tmp)
        tmp['window_start']<- ((tmp$window-1)*75000)+1
        tmp['window_end']<- tmp$window*75000
        tmp<-tmp[,c(6:9,1:5)]
        }
        }
TD_mat<-data.frame(calc_TD(slide_vcf_file_75kb@Tajima.D))
outfile_td <-paste0("data/popgenome/",tid,".td.popgenome.out")
write.table(TD_mat, outfile_td, sep = "\t", row.names=F, quote = F)

}
x_vector<-seq(1,nrow(chr_length))
lapply(x_vector,chr_function)
