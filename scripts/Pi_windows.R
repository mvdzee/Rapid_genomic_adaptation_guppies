# read in the raw data file
setwd('FIBR_ms')

pi <- read.table("data/popgenome/FIBR_Pi.txt", row.names = NULL, header = TRUE)
pi$window_start<-pi$window_start+1

# add window ID and re-order the table
pi['window_id']<-paste0(pi$chrom, ":",pi$window_start,"-",pi$window_end)
pi2<-pi[,c(1:4,16,5:15)]
pi3<-na.omit(pi2)

### create empty data frame
pi_rbind<-data.frame(matrix(ncol=7))
colnames(pi_rbind)<-c(colnames(pi3)[c(1:5)],"Pi","comp")

### loop to re-format & fill the empty dataframe. Now has a column pi and comp (=comparison)
## CHANGE 5:8 to 5:9 when UL is done!
for (i in 6:11){
  pi_tmp<-pi3[,c(1:5,i)]
  pi_tmp$comp<-rep(colnames(pi3)[i],length(pi_tmp$chrom))
  colnames(pi_tmp)[6]<-"Pi"
  pi_rbind<-rbind(pi_rbind,pi_tmp)
}
# remove NA line at the top
pi_rbind<-pi_rbind[-1,]
write.table(pi_rbind, "output/popgenome/FIBR_Pi_processed.txt", quote = F, sep = '\t', row.names = F)
