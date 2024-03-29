###########################################################
## This is an R script that processes raw Fst files created by PopGenome   
## It outputs a table with median global Fst above the diagonal and mean Fst below the diagonal  

# read in the concatenated Fst file
fst<-read.table("data/popgenome/FIBR_Fst.txt", header = T)
fst[,-c(1,2,3,4)][fst[,-c(1,2,3,4)] < 0] <- 0

# remove NA from the data columns, adapt
fst2<-na.omit(fst[,c(5:length(fst))])  

# these for loops calculate and print mean/median Fst for each pair wise population
for(i in seq(from=5,to=ncol(fst),by=1)){
  print(mean(fst[,i],na.rm = T))
  print(colnames(fst[i]))
}

for(i in seq(from=5,to=ncol(fst),by=1)){
  print(colnames(fst[i]))
  print(median(fst[,i],na.rm = T))
}

# this creates a matrix and populations it with mean (below diagonal) and median (above diagonal)
mat <- matrix(NA, nrow = 5, ncol = 5, dimnames = list(c("ILL","IUL","IC","IT","GHP"), c("ILL","IUL","IC","IT","GHP")))
#mean
mat[2,1]<-mean(fst2$ILL_IUL)
mat[3,1]<-mean(fst2$IC_ILL)
mat[4,1]<-mean(fst2$ILL_IT)
mat[5,1]<-mean(fst2$GHP_ILL)
mat[3,2]<-mean(fst2$IC_IUL)
mat[4,2]<-mean(fst2$IT_IUL)
mat[5,2]<-mean(fst2$GHP_IUL)
mat[4,3]<-mean(fst2$IC_IT)
mat[5,3]<-mean(fst2$IC_GHP)
mat[5,4]<-mean(fst2$GHP_IT)

#median
mat[1,2]<-median(fst2$ILL_IUL)
mat[1,3]<-median(fst2$IC_ILL)
mat[1,4]<-median(fst2$ILL_IT)
mat[1,5]<-median(fst2$GHP_ILL)
mat[2,3]<-median(fst2$IC_IUL)
mat[2,4]<-median(fst2$IT_IUL)
mat[2,5]<-median(fst2$GHP_IUL)
mat[3,4]<-median(fst2$IC_IT)
mat[3,5]<-median(fst2$IC_GHP)
mat[4,5]<-median(fst2$GHP_IT)

write.table(mat, "output/popgenome/fst_genomewide.txt", quote = F, sep = '\t', row.names = F)
