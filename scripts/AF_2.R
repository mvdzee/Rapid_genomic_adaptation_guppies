# this R script processes the raw AFs from vcftools. 
# it collapses the per SNP values into windows of 75kb (or size of your choice)

library(Rfast)
library(parallel)
library(data.table)

ghp_af <- fread('data/allele_freq/GHP_AF.txt', header = T)
ghp_af2<-ghp_af
ghp_af2['GHP_minor_FRQ']<-pmin(ghp_af$REF_FRQ, ghp_af$ALT_FRQ)
# Index returns 1/2 for minor allele frequency
minor_index<-rowMins(as.matrix(ghp_af[,c("REF_FRQ","ALT_FRQ")]))
ghp_af2$minor_allele<-unlist(mclapply(1:length(minor_index),function(x){return(ghp_af2[x,c("REF","ALT")[minor_index[x]]])},mc.cores = 3))
ghp_af2<-ghp_af2[with(ghp_af2,order(CHROM,POS)),]
ghp_minor<-ghp_af2[,c(1,2,9,10)]
ghp_minor$comp<-rep('GHP',nrow(ghp_minor))

### Read this in if re-doing something so you don't have to redo the function thing
write.table(ghp_minor, "output/allele_freq/GHP_minor.txt", quote = F, sep = '\t', row.names = F)
# if you nead to read in the minor alleles later
ghp_minor<-fread("output/allele_freq/GHP_minor.txt",header = T)

ll_af<-fread('data/allele_freq/ILL_AF.txt', header = T)
ll_af<-ll_af[with(ll_af, order(CHROM,POS)), ]
ll_af2<-ll_af[(ll_af$REF==ghp_minor$minor_allele),c(1,2,5,6)]
colnames(ll_af2)<-c('CHROM','POS','minor_allele','minor_freq')
ll_af3<-ll_af[(ll_af$ALT==ghp_minor$minor_allele),c(1,2,7,8)]
colnames(ll_af3)<-c('CHROM','POS','minor_allele','minor_freq')
ll_minor<-rbind(ll_af2,ll_af3)
ll_minor<-ll_minor[with(ll_minor, order(CHROM,POS)), ]
ll_minor$comp<-rep('ILL',nrow(ll_minor))
write.table(ll_minor, "output/allele_freq/LL_minor.txt", quote = F, sep = '\t', row.names = F)
# if you nead to read in the minor alleles later
ll_minor<-fread("output/allele_freq/LL_minor.txt",header=T)

ul_af<-fread('data/allele_freq/IUL_AF.txt', header = T)
ul_af<-ul_af[with(ul_af, order(CHROM,POS)), ]
ul_af2<-ul_af[(ul_af$REF==ghp_minor$minor_allele),c(1,2,5,6)]
colnames(ul_af2)<-c('CHROM','POS','minor_allele','minor_freq')
ul_af3<-ul_af[(ul_af$ALT==ghp_minor$minor_allele),c(1,2,7,8)]
colnames(ul_af3)<-c('CHROM','POS','minor_allele','minor_freq')
ul_minor<-rbind(ul_af2,ul_af3)
ul_minor<-ul_minor[with(ul_minor, order(CHROM,POS)), ]
ul_minor$comp<-rep('IUL',nrow(ul_minor))
write.table(ul_minor, "output/UL_minor.txt", quote = F, sep = '\t', row.names = F)
# if you nead to read in the minor alleles later
ul_minor<-fread("output/allele_freq/UL_minor.txt",header=T)

c_af<-fread('data/allele_freq/IC_AF.txt', header = T)
c_af<-c_af[with(c_af, order(CHROM,POS)), ]
c_af2<-c_af[(c_af$REF==ghp_minor$minor_allele),c(1,2,5,6)]
colnames(c_af2)<-c('CHROM','POS','minor_allele','minor_freq')
c_af3<-c_af[(c_af$ALT==ghp_minor$minor_allele),c(1,2,7,8)]
colnames(c_af3)<-c('CHROM','POS','minor_allele','minor_freq')
c_minor<-rbind(c_af2,c_af3)
c_minor<-c_minor[with(c_minor, order(CHROM,POS)), ]
c_minor$comp<-rep('IC',nrow(c_minor))
write.table(c_minor, "output/allele_freq/C_minor.txt", quote = F, sep = '\t', row.names = F)
# if you nead to read in the minor alleles later
c_minor<-fread("output/allele_freq/C_minor.txt",header=T)

t_af<-fread('data/allele_freq/IT_AF.txt', header = T)
t_af<-t_af[with(t_af, order(CHROM,POS)), ]
t_af2<-t_af[(t_af$REF==ghp_minor$minor_allele),c(1,2,5,6)]
colnames(t_af2)<-c('CHROM','POS','minor_allele','minor_freq')
t_af3<-t_af[(t_af$ALT==ghp_minor$minor_allele),c(1,2,7,8)]
colnames(t_af3)<-c('CHROM','POS','minor_allele','minor_freq')
t_minor<-rbind(t_af2,t_af3)
t_minor<-t_minor[with(t_minor, order(CHROM,POS)), ]
t_minor$comp<-rep('IT',nrow(t_minor))
write.table(t_minor, "output/allele_freq/T_minor.txt", quote = F, sep = '\t', row.names = F)
# if you nead to read in the minor alleles later
t_minor<-fread("output/allele_freq/T_minor.txt",header=T)


### calculate delta AFs
# intro pairs
dLL <- data.frame(ll_minor$minor_freq - ghp_minor$GHP_minor_FRQ)
dLL['chrom']<-ll_minor$CHROM
dLL['BP']<-ll_minor$POS
dLL['comp']<-rep('ILL',nrow(dLL))
dLL<-dLL[,c(2,3,1,4)]
colnames(dLL)<-c('chrom','BP','delta_AF','comp')
dLL$abs_AF<-abs(dLL$delta_AF)

dUL <- data.frame(ul_minor$minor_freq - ghp_minor$GHP_minor_FRQ)
dUL['chrom']<-ul_minor$CHROM
dUL['BP']<-ul_minor$POS
dUL['comp']<-rep('IUL',nrow(dUL))
dUL<-dUL[,c(2,3,1,4)]
colnames(dUL)<-c('chrom','BP','delta_AF','comp')
dUL$abs_AF<-abs(dUL$delta_AF)

dC <- data.frame(c_minor$minor_freq - ghp_minor$GHP_minor_FRQ)
dC['chrom']<-c_minor$CHROM
dC['BP']<-c_minor$POS
dC['comp']<-rep('IC',nrow(dC))
dC<-dC[,c(2,3,1,4)]
colnames(dC)<-c('chrom','BP','delta_AF','comp')
dC$abs_AF<-abs(dC$delta_AF)

dT <- data.frame(t_minor$minor_freq - ghp_minor$GHP_minor_FRQ)
dT['chrom']<-t_minor$CHROM
dT['BP']<-t_minor$POS
dT['comp']<-rep('IT',nrow(dT))
dT<-dT[,c(2,3,1,4)]
colnames(dT)<-c('chrom','BP','delta_AF','comp')
dT$abs_AF<-abs(dT$delta_AF)

AF_all<-rbind(dLL,dUL,dC,dT)
write.table(AF_all, "output/allele_freq/DAF_perSNP.txt", quote = F, sep = '\t', row.names = F)

# Average into windows non-absolute windows
chrs<-unique(AF_all$chrom)
wind_size<-75000

### This will average the delta AFs into windows
winds_LL<-data.frame(rbindlist(mclapply(chrs,function(x){
  tmp<-dLL[dLL$chrom == x,]
  winds1<-seq(0,max(tmp$BP),by=wind_size)
  winds2<-winds1+wind_size
  
  # Summarise for each
  sum_AF<-data.frame(rbindlist(lapply(1:length(winds2),function(y){
    tmp2<-tmp[tmp$BP <= winds2[y] & tmp$BP >= winds1[y],]
    out<-data.frame(mean(tmp2$delta_AF,na.rm = T))
    # Tidy
    #out$river<-rownames(out)
    out$chrom<-x
    out$window<-y
    out$BP1<-as.integer(winds1[y])+1
    out$BP2<-as.integer(winds2[y])
    out$comp<-'ILL'
    colnames(out)<-c('delta_AF','chrom','window','BP1','BP2','comp')
    
    return(out)
  })))
  return(sum_AF)
},mc.cores=3)))


winds_UL<-data.frame(rbindlist(mclapply(chrs,function(x){
  tmp<-dUL[dUL$chrom == x,]
  winds1<-seq(0,max(tmp$BP),by=wind_size)
  winds2<-winds1+wind_size
  
  # Summarise for each
  sum_AF<-data.frame(rbindlist(lapply(1:length(winds2),function(y){
    tmp2<-tmp[tmp$BP <= winds2[y] & tmp$BP >= winds1[y],]
    #out<-data.frame(dLL_chr=colMeans(abs(dLL_chr[,3])))
    out<-data.frame(mean(tmp2$delta_AF))
    # Tidy
    #out$river<-rownames(out)
    out$chrom<-x
    out$window<-y
    out$BP1<-as.integer(winds1[y])+1
    out$BP2<-as.integer(winds2[y])
    out$comp<-'IUL'
    colnames(out)<-c('delta_AF','chrom','window','BP1','BP2','comp')
    return(out)
  })))
  return(sum_AF)
},mc.cores=3)))

winds_C<-data.frame(rbindlist(mclapply(chrs,function(x){
  tmp<-dC[dC$chrom == x,]
  winds1<-seq(0,max(tmp$BP),by=wind_size)
  winds2<-winds1+wind_size
  
  # Summarise for each
  sum_AF<-data.frame(rbindlist(lapply(1:length(winds2),function(y){
    tmp2<-tmp[tmp$BP <= winds2[y] & tmp$BP >= winds1[y],]
    #out<-data.frame(dLL_chr=colMeans(abs(dLL_chr[,3])))
    out<-data.frame(mean(tmp2$delta_AF))
    # Tidy
    #out$river<-rownames(out)
    out$chrom<-x
    out$window<-y
    out$BP1<-as.integer(winds1[y])+1
    out$BP2<-as.integer(winds2[y])
    out$comp<-'IC'
    colnames(out)<-c('delta_AF','chrom','window','BP1','BP2','comp')
    return(out)
  })))
  return(sum_AF)
},mc.cores=3)))

winds_T<-data.frame(rbindlist(mclapply(chrs,function(x){
  tmp<-dT[dT$chrom == x,]
  winds1<-seq(0,max(tmp$BP),by=wind_size)
  winds2<-winds1+wind_size
  
  # Summarise for each
  sum_AF<-data.frame(rbindlist(lapply(1:length(winds2),function(y){
    tmp2<-tmp[tmp$BP <= winds2[y] & tmp$BP >= winds1[y],]
    #out<-data.frame(dLL_chr=colMeans(abs(dLL_chr[,3])))
    out<-data.frame(mean(tmp2$delta_AF))
    # Tidy
    #out$river<-rownames(out)
    out$chrom<-x
    out$window<-y
    out$BP1<-as.integer(winds1[y])+1
    out$BP2<-as.integer(winds2[y])
    out$comp<-'IT'
    colnames(out)<-c('delta_AF','chrom','window','BP1','BP2','comp')
    return(out)
  })))
  return(sum_AF)
},mc.cores=3)))

delta_AF <-rbind(winds_LL,winds_UL,winds_C,winds_T)
write.table(delta_AF, "output/allele_freq/DAF_NonAbs_windows.txt", quote = F, sep = '\t', row.names = F)


### This will average absolute delta AFs into windows
winds_LL2<-data.frame(rbindlist(mclapply(chrs,function(x){
  tmp<-dLL[dLL$chrom == x,]
  winds1<-seq(0,max(tmp$BP),by=wind_size)
  winds2<-winds1+wind_size
  
  # Summarise for each
  sum_AF<-data.frame(rbindlist(lapply(1:length(winds2),function(y){
    tmp2<-tmp[tmp$BP <= winds2[y] & tmp$BP >= winds1[y],]
    out<-data.frame(mean(tmp2$abs_AF,na.rm = T))
    # Tidy
    #out$river<-rownames(out)
    out$chrom<-x
    out$window<-y
    out$BP1<-as.integer(winds1[y])+1
    out$BP2<-as.integer(winds2[y])
    out$comp<-'ILL'
    colnames(out)<-c('abs_AF','chrom','window','BP1','BP2','comp')
    
    return(out)
  })))
  return(sum_AF)
},mc.cores=3)))


winds_UL2<-data.frame(rbindlist(mclapply(chrs,function(x){
  tmp<-dUL[dUL$chrom == x,]
  winds1<-seq(0,max(tmp$BP),by=wind_size)
  winds2<-winds1+wind_size
  
  # Summarise for each
  sum_AF<-data.frame(rbindlist(lapply(1:length(winds2),function(y){
    tmp2<-tmp[tmp$BP <= winds2[y] & tmp$BP >= winds1[y],]
    #out<-data.frame(dLL_chr=colMeans(abs(dLL_chr[,3])))
    out<-data.frame(mean(tmp2$abs_AF))
    # Tidy
    #out$river<-rownames(out)
    out$chrom<-x
    out$window<-y
    out$BP1<-as.integer(winds1[y])+1
    out$BP2<-as.integer(winds2[y])
    out$comp<-'IUL'
    colnames(out)<-c('abs_AF','chrom','window','BP1','BP2','comp')
    return(out)
  })))
  return(sum_AF)
},mc.cores=3)))

winds_C2<-data.frame(rbindlist(mclapply(chrs,function(x){
  tmp<-dC[dC$chrom == x,]
  winds1<-seq(0,max(tmp$BP),by=wind_size)
  winds2<-winds1+wind_size
  
  # Summarise for each
  sum_AF<-data.frame(rbindlist(lapply(1:length(winds2),function(y){
    tmp2<-tmp[tmp$BP <= winds2[y] & tmp$BP >= winds1[y],]
    #out<-data.frame(dLL_chr=colMeans(abs(dLL_chr[,3])))
    out<-data.frame(mean(tmp2$abs_AF))
    # Tidy
    #out$river<-rownames(out)
    out$chrom<-x
    out$window<-y
    out$BP1<-as.integer(winds1[y])+1
    out$BP2<-as.integer(winds2[y])
    out$comp<-'IC'
    colnames(out)<-c('abs_AF','chrom','window','BP1','BP2','comp')
    return(out)
  })))
  return(sum_AF)
},mc.cores=3)))

winds_T2<-data.frame(rbindlist(mclapply(chrs,function(x){
  tmp<-dT[dT$chrom == x,]
  winds1<-seq(0,max(tmp$BP),by=wind_size)
  winds2<-winds1+wind_size
  
  # Summarise for each
  sum_AF<-data.frame(rbindlist(lapply(1:length(winds2),function(y){
    tmp2<-tmp[tmp$BP <= winds2[y] & tmp$BP >= winds1[y],]
    #out<-data.frame(dLL_chr=colMeans(abs(dLL_chr[,3])))
    out<-data.frame(mean(tmp2$abs_AF))
    # Tidy
    #out$river<-rownames(out)
    out$chrom<-x
    out$window<-y
    out$BP1<-as.integer(winds1[y])+1
    out$BP2<-as.integer(winds2[y])
    out$comp<-'IT'
    colnames(out)<-c('abs_AF','chrom','window','BP1','BP2','comp')
    return(out)
  })))
  return(sum_AF)
},mc.cores=3)))

delta_AF2 <-rbind(winds_LL2,winds_UL2,winds_C2,winds_T2)
write.table(delta_AF2, "output/allele_freq/DAF_Abs_windows.txt", quote = F, sep = '\t', row.names = F)
