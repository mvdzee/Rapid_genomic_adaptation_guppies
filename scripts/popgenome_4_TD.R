library(data.table)

### This script processes the raw Tajimas D values from popgenome
td <- read.table("data/td/FIBR_TD.txt", row.names = NULL, header = TRUE)
td['window_id']<-paste0(td$chrom, ":",td$window_start,"-",td$window_end)
td2<-na.omit(td2)

### create empty data frame
td_rbind<-data.frame(matrix(ncol=7))
colnames(td_rbind)<-c(colnames(td2)[c(1:5)],"TD","comp")

### loop to re-format & fill the empty dataframe. Now has a column pi and comp (=comparison)
## CHANGE 5:8 to 5:9 when UL is done!
for (i in 6:10){
  td_tmp<-td2[,c(1:5,i)]
  td_tmp$comp<-rep(colnames(td2)[i],length(td_tmp$chrom))
  colnames(td_tmp)[6]<-"TD"
  td_rbind<-rbind(td_rbind,td_tmp)
}
# remove NA line at the top
td_rbind<-td_rbind[-1,]
write.table(td_rbind, "output/td/FIBR_TD_processed.txt", quote = F, sep = '\t', row.names = F)

