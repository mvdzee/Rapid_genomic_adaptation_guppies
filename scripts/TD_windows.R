# This R script processes the raw, concatenated Tajima's D file created by PopGenome

library(data.table)
library(stringr)
library(VennDiagram)

# Read in the concatenated file and reorder the table
td <- read.table("data/popgenome/FIBR_TD.txt", row.names = NULL, header = TRUE)
td$window_start<-td$window_start+1
td['window_id']<-paste0(td$chrom, ":",td$window_start,"-",td$window_end)
td2<-td[,c(1:4,11,5:10)]
td2<-na.omit(td2)

### create empty data frame
td_rbind<-data.frame(matrix(ncol=7))
colnames(td_rbind)<-c(colnames(td2)[c(1:5)],"TD","comp")

### loop to re-format & fill the empty dataframe. Now has a column pi and comp (=comparison)
## Change '6:11' to include the appropriate columns from your data
for (i in 6:11){
  pi_tmp<-td2[,c(1:5,i)]
  pi_tmp$comp<-rep(colnames(td2)[i],length(pi_tmp$chrom))
  colnames(pi_tmp)[6]<-"TD"
  td_rbind<-rbind(td_rbind,pi_tmp)
}
# remove NA line at the top
td_rbind<-td_rbind[-1,]
write.table(td_rbind, "output/popgenome/FIBR_TD_processed.txt", quote = F, sep = '\t', row.names = F)
