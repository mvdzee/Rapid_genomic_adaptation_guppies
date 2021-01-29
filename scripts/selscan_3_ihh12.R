library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(parallel)

### read in all gle chromosomes
c<-fread("data/selscan/ihh12/C_ihh12.norm.txt",header = T,fill=TRUE)
t<-fread("data/selscan/ihh12/T_ihh12.norm.txt",header = T,fill=TRUE)
ll<-fread("data/selscan/ihh12/LL_ihh12.norm.txt",header = T,fill=TRUE)
ul<-fread("data/selscan/ihh12/UL_ihh12.norm.txt",header = T,fill=TRUE)
gh<-fread("data/selscan/ihh12/GH_ihh12.norm.txt",header = T,fill=TRUE)

c$POP<-rep("IC",nrow(c))
t$POP<-rep("IT",nrow(t))
ll$POP<-rep("ILL",nrow(ll))
ul$POP<-rep("IUL",nrow(ul))
gh$POP<-rep("GHP",nrow(gh))

c<-na.omit(c)
t<-na.omit(t)
ll<-na.omit(ll)
ul<-na.omit(ul)
gh<-na.omit(gh)

c<-separate(data=c,col=id, into = c("chrom","pos"),sep="\\:",remove = F)
t<-separate(data=t,col=id, into = c("chrom","pos"),sep="\\:",remove = F)
ll<-separate(data=ll,col=id, into = c("chrom","pos"),sep="\\:",remove = F)
ul<-separate(data=ul,col=id, into = c("chrom","pos"),sep="\\:",remove = F)
gh<-separate(data=gh,col=id, into = c("chrom","pos"),sep="\\:",remove = F)

c$pos<-as.numeric(c$pos)
t$pos<-as.numeric(t$pos)
ll$pos<-as.numeric(ll$pos)
ul$pos<-as.numeric(ul$pos)
gh$pos<-as.numeric(gh$pos)

# average into windows
chrs<-unique(c$chrom)
wind_size<-75000
winds_c<-data.frame(rbindlist(mclapply(chrs,function(x){
  tmp<-c[c$chrom == x,]
  if(nrow(tmp) == 0){
    print("empty")
  }else{
    winds1<-seq(0,max(tmp$pos),by=wind_size)
    winds2<-winds1+wind_size
  }
  # Summarise for each
  sum_fst<-data.frame(rbindlist(mclapply(1:length(winds2),function(y){
    tmp2<-tmp[which(tmp$pos>=winds1[y] & tmp$pos<=winds2[y]),]
    if(nrow(tmp2) == 0){
      out <- data.frame("NA")
    } else {
      out<-data.frame(as.numeric(mean(tmp2$normihh12,na.rm = T)))
    }
    # Tidy
    out$chrom<-x
    out$window<-y
    out$BP1<-as.integer(winds1[y])+1
    out$BP2<-as.integer(winds2[y])
    out$comp<-'IC'
    colnames(out)<-c('normihh12','chrom',"window",'BP1','BP2','comp')
    out<-out[,c(2:5,1,6)]
    return(out)
  },mc.cores = 3)))
  return(sum_fst)
},mc.cores = 3)))
winds_c$normihh12<-as.numeric(as.character(winds_c$normihh12))
winds_c$ihh12<-as.numeric(as.character(winds_c$ihh12))
winds_c$window_id<-paste0(winds_c$chrom,':',winds_c$BP1,'-',winds_c$BP2)
write.table(winds_c, "output/selscan/ihh12/IC_norm.ihh12_windows.txt", quote = F, sep = '\t', row.names = F)

chrs<-unique(t$chrom)
wind_size<-75000
winds_t<-data.frame(rbindlist(mclapply(chrs,function(x){
  tmp<-t[t$chrom == x,]
  if(nrow(tmp) == 0){
    print("empty")
  }else{
    winds1<-seq(0,max(tmp$pos),by=wind_size)
    winds2<-winds1+wind_size
  }
  # Summarise for each
  sum_fst<-data.frame(rbindlist(mclapply(1:length(winds2),function(y){
    tmp2<-tmp[which(tmp$pos>=winds1[y] & tmp$pos<=winds2[y]),]
    if(nrow(tmp2) == 0){
      out <- data.frame("NA")
    } else {
      out<-data.frame(as.numeric(mean(tmp2$normihh12,na.rm = T)))
    }
    # Tidy
    out$chrom<-x
    out$window<-y
    out$BP1<-as.integer(winds1[y])+1
    out$BP2<-as.integer(winds2[y])
    out$comp<-'IT'
    colnames(out)<-c('normihh12','chrom',"window",'BP1','BP2','comp')
    out<-out[,c(2:5,1,6)]
    return(out)
  },mc.cores = 3)))
  return(sum_fst)
},mc.cores = 3)))
winds_t$normihh12<-as.numeric(as.character(winds_t$normihh12))
winds_t$ihh12<-as.numeric(as.character(winds_t$ihh12))
winds_t$window_id<-paste0(winds_t$chrom,':',winds_t$BP1,'-',winds_t$BP2)
write.table(winds_t, "output/selscan/ihh12/IT_norm.ihh12_windows.txt", quote = F, sep = '\t', row.names = F)

chrs<-unique(ll$chrom)
wind_size<-75000
winds_ll<-data.frame(rbindlist(mclapply(chrs,function(x){
  tmp<-ll[ll$chrom == x,]
  if(nrow(tmp) == 0){
    print("empty")
  }else{
    winds1<-seq(0,max(tmp$pos),by=wind_size)
    winds2<-winds1+wind_size
  }
  # Summarise for each
  sum_fst<-data.frame(rbindlist(mclapply(1:length(winds2),function(y){
    tmp2<-tmp[which(tmp$pos>=winds1[y] & tmp$pos<=winds2[y]),]
    if(nrow(tmp2) == 0){
      out <- data.frame("NA")
    } else {
      out<-data.frame(as.numeric(mean(tmp2$normihh12,na.rm = T)))
    }
    
    # Tidy
    out$chrom<-x
    out$window<-y
    out$BP1<-as.integer(winds1[y])+1
    out$BP2<-as.integer(winds2[y])
    out$comp<-'ILL'
    colnames(out)<-c('normihh12','chrom',"window",'BP1','BP2','comp')
    out<-out[,c(2:5,1,6)]
    return(out)
  },mc.cores = 3)))
  return(sum_fst)
},mc.cores = 3)))
winds_ll$normihh12<-as.numeric(as.character(winds_ll$normihh12))
winds_ll$ihh12<-as.numeric(as.character(winds_ll$ihh12))
winds_ll$window_id<-paste0(winds_ll$chrom,':',winds_ll$BP1,'-',winds_ll$BP2)
write.table(winds_ll, "output/selscan/ihh12/ILL_norm.ihh12_windows.txt", quote = F, sep = '\t', row.names = F)

chrs<-unique(ul$chrom)
wind_size<-75000
winds_ul<-data.frame(rbindlist(mclapply(chrs,function(x){
  tmp<-ul[ul$chrom == x,]
  if(nrow(tmp) == 0){
    print("empty")
  }else{
    winds1<-seq(0,max(tmp$pos),by=wind_size)
    winds2<-winds1+wind_size
  }
  # Summarise for each
  sum_fst<-data.frame(rbindlist(mclapply(1:length(winds2),function(y){
    tmp2<-tmp[which(tmp$pos>=winds1[y] & tmp$pos<=winds2[y]),]
    if(nrow(tmp2) == 0){
      out <- data.frame("NA")
    } else {
      out<-data.frame(as.numeric(mean(tmp2$normihh12,na.rm = T)))
    }
    
    # Tidy
    out$chrom<-x
    out$window<-y
    out$BP1<-as.integer(winds1[y])+1
    out$BP2<-as.integer(winds2[y])
    out$comp<-'IUL'
    colnames(out)<-c('normihh12','chrom',"window",'BP1','BP2','comp')
    out<-out[,c(2:5,1,6)]
    return(out)
  },mc.cores = 3)))
  return(sum_fst)
},mc.cores = 3)))
winds_ul$normihh12<-as.numeric(as.character(winds_ul$normihh12))
winds_ul$ihh12<-as.numeric(as.character(winds_ul$ihh12))
winds_ul$window_id<-paste0(winds_ul$chrom,':',winds_ul$BP1,'-',winds_ul$BP2)
write.table(winds_ul, "output/selscan/ihh12/IUL_norm.ihh12_windows.txt", quote = F, sep = '\t', row.names = F)


chrs<-unique(gh$chrom)
wind_size<-75000
winds_gh<-data.frame(rbindlist(mclapply(chrs,function(x){
  tmp<-gh[gh$chrom == x,]
  if(nrow(tmp) == 0){
    print("empty")
  }else{
    winds1<-seq(0,max(tmp$pos),by=wind_size)
    winds2<-winds1+wind_size
  }
  # Summarise for each
  sum_fst<-data.frame(rbindlist(mclapply(1:length(winds2),function(y){
    tmp2<-tmp[which(tmp$pos>=winds1[y] & tmp$pos<=winds2[y]),]
    if(nrow(tmp2) == 0){
      out <- data.frame("NA")
    } else {
      out<-data.frame(as.numeric(mean(tmp2$normihh12,na.rm = T)))
    }
    
    # Tidy
    out$chrom<-x
    out$window<-y
    out$BP1<-as.integer(winds1[y])+1
    out$BP2<-as.integer(winds2[y])
    out$comp<-'GHP'
    colnames(out)<-c('normihh12','chrom',"window",'BP1','BP2','comp')
    out<-out[,c(2:5,1,6)]
    return(out)
  },mc.cores = 3)))
  return(sum_fst)
},mc.cores = 3)))
winds_gh$normihh12<-as.numeric(as.character(winds_gh$normihh12))
winds_gh$window_id<-paste0(winds_gh$chrom,':',winds_gh$BP1,'-',winds_gh$BP2)
write.table(winds_gh, "output/selscan/ihh12/GHP_norm.ihh12_windows.txt", quote = F, sep = '\t', row.names = F)

### process windowed data
# read in the windowed data
gh2<-na.omit(read.table('output/selscan/ihh12/GHP_norm.ihh12_windows.txt',header = T))
c2<-na.omit(read.table('output/selscan/ihh12/IC_norm.ihh12_windows.txt',header = T))
t2<-na.omit(read.table('output/selscan/ihh12/IT_norm.ihh12_windows.txt',header = T))
ll2<-na.omit(read.table('output/selscan/ihh12/ILL_norm.ihh12_windows.txt',header = T))
ul2<-na.omit(read.table('output/selscan/ihh12/IUL_norm.ihh12_windows.txt',header = T))
all<-rbind(gh2,c2,t2,ll2,ul2)
write.table(all, "output/selscan/ihh12/FIBR_ihh12_all.txt", quote = F, sep = '\t', row.names = F)

gh2$absihh12<-abs(gh2$normihh12)
c2$absihh12<-abs(c2$normihh12)
t2$absihh12<-abs(t2$normihh12)
ll2$absihh12<-abs(ll2$normihh12)
ul2$absihh12<-abs(ul2$normihh12)

# find the outliers
# create a df for each intro and GHP
c3<-c2[c2$window_id%in%gh2$window_id,]
gh3<-gh2[gh2$window_id%in%c3$window_id,]
t3<-t2[t2$window_id%in%gh2$window_id,]
gh3t<-gh2[gh2$window_id%in%t3$window_id,]
ll3<-ll2[ll2$window_id%in%gh2$window_id,]
gh3ll<-gh2[gh2$window_id%in%ll3$window_id,]
ul3<-ul2[ul2$window_id%in%gh2$window_id,]
gh3ul<-gh2[gh2$window_id%in%ul3$window_id,]

# find the outliers that are >5 in intro and <5 in GHP
ghp_c<-cbind(c3,gh3$absihh12)
ghp_c$comp2<-'IC_GHP'
ghp_c<-ghp_c[,c(1:4,8:11)]
colnames(ghp_c)<-c('chrom','window','BP1','BP2','ID','LP_ihh','GH_ihh','comp')
ghp_c_out<-ghp_c[ghp_c$LP_ihh>5 & ghp_c$GH_ihh<5,]
write.table(ghp_c_out, "output/ihh12/C_GH_outliers5.txt", quote = F, sep = '\t', row.names = F)

ghp_t<-cbind(t3,gh3t$absihh12)
ghp_t$comp2<-'IT_GHP'
ghp_t<-ghp_t[,c(1:4,8:11)]
colnames(ghp_t)<-c('chrom','window','BP1','BP2','ID','LP_ihh','GH_ihh','comp')
ghp_t_out<-ghp_t[ghp_t$LP_ihh>5 & ghp_t$GH_ihh<5,]
write.table(ghp_t_out, "output/ihh12/T_GH_outliers5.txt", quote = F, sep = '\t', row.names = F)

ghp_ll<-cbind(ll3,gh3ll$absihh12)
ghp_ll$comp2<-'ILL_GHP'
ghp_ll<-ghp_ll[,c(1:4,8:11)]
colnames(ghp_ll)<-c('chrom','window','BP1','BP2','ID','LP_ihh','GH_ihh','comp')
ghp_ll_out<-ghp_ll[ghp_ll$LP_ihh>5 & ghp_ll$GH_ihh<5,]
write.table(ghp_ll_out, "output/ihh12/LL_GH_outliers5.txt", quote = F, sep = '\t', row.names = F)

ghp_ul<-cbind(ul3,gh3ul$absihh12)
ghp_ul$comp2<-'IUL_GHP'
ghp_ul<-ghp_ul[,c(1:4,8:11)]
colnames(ghp_ul)<-c('chrom','window','BP1','BP2','ID','LP_ihh','GH_ihh','comp')
ghp_ul_out<-ghp_ul[ghp_ul$LP_ihh>5 & ghp_ul$GH_ihh<5,]
write.table(ghp_ul_out, "output/ihh12/UL_GH_outliers5.txt", quote = F, sep = '\t', row.names = F)

comp_all<-rbind(ghp_c_out,ghp_t_out,ghp_ll_out,ghp_ul_out)

# save all outliers
write.table(comp_all, "output/selscan/ihh12/FIBR_ihh12_GH_LP_outliers.txt", quote = F, sep = '\t', row.names = F)

