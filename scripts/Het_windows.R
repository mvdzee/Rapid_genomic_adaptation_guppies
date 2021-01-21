library(data.table)
library(parallel)
setwd('FIBR_ms')

c_het<-read.table("data/heterozygosity/IC_het.txt",header = T)
c_het<-c_het[,c(1,2,6,7,8)]
colnames(c_het)<-c("CHR","POS","HOM1","HET","HOM2")
c_het$exp_het<-(c_het$HET/(c_het$HOM1+c_het$HET+c_het$HOM2))
c_het<-c_het[,c(1,2,6)]

t_het<-read.table("data/heterozygosity/IT_het.txt",header = T)
t_het<-t_het[,c(1,2,6,7,8)]
colnames(t_het)<-c("CHR","POS","HOM1","HET","HOM2")
t_het$exp_het<-(t_het$HET/(t_het$HOM1+t_het$HET+t_het$HOM2))
t_het<-t_het[,c(1,2,6)]

ll_het<-read.table("data/heterozygosity/ILL_het.txt",header = T)
ll_het<-ll_het[,c(1,2,6,7,8)]
colnames(ll_het)<-c("CHR","POS","HOM1","HET","HOM2")
ll_het$exp_het<-(ll_het$HET/(ll_het$HOM1+ll_het$HET+ll_het$HOM2))
ll_het<-ll_het[,c(1,2,6)]

ul_het<-read.table("data/heterozygosity/IUL_het.txt",header = T)
ul_het<-ul_het[,c(1,2,6,7,8)]
colnames(ul_het)<-c("CHR","POS","HOM1","HET","HOM2")
ul_het$exp_het<-(ul_het$HET/(ul_het$HOM1+ul_het$HET+ul_het$HOM2))
ul_het<-ul_het[,c(1,2,6)]

ghp_het<-read.table("data/heterozygosity/GHP_het.txt",header = T)
ghp_het<-ghp_het[,c(1,2,6,7,8)]
colnames(ghp_het)<-c("CHR","POS","HOM1","HET","HOM2")
ghp_het$exp_het<-(ghp_het$HET/(ghp_het$HOM1+ghp_het$HET+ghp_het$HOM2))
ghp_het<-ghp_het[,c(1,2,6)]

glp_het<-read.table("data/heterozygosity/GLP_het.txt",header = T)
glp_het<-glp_het[,c(1,2,6,7,8)]
colnames(glp_het)<-c("CHR","POS","HOM1","HET","HOM2")
glp_het$exp_het<-(glp_het$HET/(glp_het$HOM1+glp_het$HET+glp_het$HOM2))
glp_het<-glp_het[,c(1,2,6)]

## average the per SNP values into windows
chrs<-unique(c_het$CHR)
wind_size<-75000

winds_c<-data.frame(rbindlist(mclapply(chrs,function(x){
  tmp<-c_het[c_het$CHR == x,]
  winds1<-seq(0,max(tmp$POS),by=wind_size)
  winds2<-winds1+wind_size
  
  # Summarise for each
  sum_AF<-data.frame(rbindlist(lapply(1:length(winds2),function(y){
    tmp2<-tmp[tmp$POS <= winds2[y] & tmp$POS >= winds1[y],]
    out<-data.frame(mean(tmp2$exp_het,na.rm = T))
    # Tidy
    #out$river<-rownames(out)
    out$chrom<-x
    out$window<-y
    out$BP1<-as.integer(winds1[y])+1
    out$BP2<-as.integer(winds2[y])
    out$comp<-'IC'
    colnames(out)<-c('exp_het','chrom','window','BP1','BP2','comp')
    
    return(out)
  })))
  return(sum_AF)
},mc.cores=3)))

winds_t<-data.frame(rbindlist(mclapply(chrs,function(x){
  tmp<-t_het[t_het$CHR == x,]
  winds1<-seq(0,max(tmp$POS),by=wind_size)
  winds2<-winds1+wind_size
  
  # Summarise for each
  sum_AF<-data.frame(rbindlist(lapply(1:length(winds2),function(y){
    tmp2<-tmp[tmp$POS <= winds2[y] & tmp$POS >= winds1[y],]
    out<-data.frame(mean(tmp2$exp_het,na.rm = T))
    # Tidy
    #out$river<-rownames(out)
    out$chrom<-x
    out$window<-y
    out$BP1<-as.integer(winds1[y])+1
    out$BP2<-as.integer(winds2[y])
    out$comp<-'IT'
    colnames(out)<-c('exp_het','chrom','window','BP1','BP2','comp')
    
    return(out)
  })))
  return(sum_AF)
},mc.cores=3)))

winds_ll<-data.frame(rbindlist(mclapply(chrs,function(x){
  tmp<-ll_het[ll_het$CHR == x,]
  winds1<-seq(0,max(tmp$POS),by=wind_size)
  winds2<-winds1+wind_size
  
  # Summarise for each
  sum_AF<-data.frame(rbindlist(lapply(1:length(winds2),function(y){
    tmp2<-tmp[tmp$POS <= winds2[y] & tmp$POS >= winds1[y],]
    out<-data.frame(mean(tmp2$exp_het,na.rm = T))
    # Tidy
    #out$river<-rownames(out)
    out$chrom<-x
    out$window<-y
    out$BP1<-as.integer(winds1[y])+1
    out$BP2<-as.integer(winds2[y])
    out$comp<-'ILL'
    colnames(out)<-c('exp_het','chrom','window','BP1','BP2','comp')
    
    return(out)
  })))
  return(sum_AF)
},mc.cores=3)))


winds_ul<-data.frame(rbindlist(mclapply(chrs,function(x){
  tmp<-ul_het[ul_het$CHR == x,]
  winds1<-seq(0,max(tmp$POS),by=wind_size)
  winds2<-winds1+wind_size
  
  # Summarise for each
  sum_AF<-data.frame(rbindlist(lapply(1:length(winds2),function(y){
    tmp2<-tmp[tmp$POS <= winds2[y] & tmp$POS >= winds1[y],]
    out<-data.frame(mean(tmp2$exp_het,na.rm = T))
    # Tidy
    #out$river<-rownames(out)
    out$chrom<-x
    out$window<-y
    out$BP1<-as.integer(winds1[y])+1
    out$BP2<-as.integer(winds2[y])
    out$comp<-'IUL'
    colnames(out)<-c('exp_het','chrom','window','BP1','BP2','comp')
    
    return(out)
  })))
  return(sum_AF)
},mc.cores=3)))


winds_ghp<-data.frame(rbindlist(mclapply(chrs,function(x){
  tmp<-ghp_het[ghp_het$CHR == x,]
  winds1<-seq(0,max(tmp$POS),by=wind_size)
  winds2<-winds1+wind_size
  
  # Summarise for each
  sum_AF<-data.frame(rbindlist(lapply(1:length(winds2),function(y){
    tmp2<-tmp[tmp$POS <= winds2[y] & tmp$POS >= winds1[y],]
    out<-data.frame(mean(tmp2$exp_het,na.rm = T))
    # Tidy
    #out$river<-rownames(out)
    out$chrom<-x
    out$window<-y
    out$BP1<-as.integer(winds1[y])+1
    out$BP2<-as.integer(winds2[y])
    out$comp<-'GHP'
    colnames(out)<-c('exp_het','chrom','window','BP1','BP2','comp')
    
    return(out)
  })))
  return(sum_AF)
},mc.cores=3)))

winds_glp<-data.frame(rbindlist(mclapply(chrs,function(x){
  tmp<-glp_het[glp_het$CHR == x,]
  winds1<-seq(0,max(tmp$POS),by=wind_size)
  winds2<-winds1+wind_size
  
  # Summarise for each
  sum_AF<-data.frame(rbindlist(lapply(1:length(winds2),function(y){
    tmp2<-tmp[tmp$POS <= winds2[y] & tmp$POS >= winds1[y],]
    out<-data.frame(mean(tmp2$exp_het,na.rm = T))
    # Tidy
    #out$river<-rownames(out)
    out$chrom<-x
    out$window<-y
    out$BP1<-as.integer(winds1[y])+1
    out$BP2<-as.integer(winds2[y])
    out$comp<-'GLP'
    colnames(out)<-c('exp_het','chrom','window','BP1','BP2','comp')
    
    return(out)
  })))
  return(sum_AF)
},mc.cores=3)))


het_all <-rbind(winds_ghp,winds_glp,winds_ll,winds_ul,winds_c,winds_t)
write.table(het_all, "output/heterozygosity/FIBR_heterozygosity.txt", quote = F, sep = '\t', row.names = F)

