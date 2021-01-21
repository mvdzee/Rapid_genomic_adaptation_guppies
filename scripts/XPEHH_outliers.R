library(tidyr)
library(data.table)
library(ggplot2)
library(parallel)
library(cowplot)

# read in the concatenated files
ghp_c<-fread("data/selscan/xpehh/GH_C_XPEHH.txt",header = T)
ghp_t<-fread("data/selscan/xpehh/GH_T_XPEHH.txt",header = T)
ghp_ll<-fread("data/selscan/xpehh/GH_LL_XPEHH.txt",header = T)
ghp_ul<-fread("data/selscan/xpehh/GH_UL_XPEHH.txt",header = T)

#ghp_glp<-separate(data=ghp_glp,col=id, into = c("chrom","pos"),sep="\\:",remove=F)
ghp_c<-separate(data=ghp_c,col=id, into = c("chrom","pos"),sep="\\:",remove=F)
ghp_t<-separate(data=ghp_t,col=id, into = c("chrom","pos"),sep="\\:",remove=F)
ghp_ll<-separate(data=ghp_ll,col=id, into = c("chrom","pos"),sep="\\:",remove=F)
ghp_ul<-separate(data=ghp_ul,col=id, into = c("chrom","pos"),sep="\\:",remove=F)

#ghp_glp$pos<-as.numeric(ghp_glp$pos)
ghp_c$pos<-as.numeric(ghp_c$pos)
ghp_t$pos<-as.numeric(ghp_t$pos)
ghp_ll$pos<-as.numeric(ghp_ll$pos)
ghp_ul$pos<-as.numeric(ghp_ul$pos)

# collapse into windows
chrs<-unique(ghp_c$chrom)
wind_size<-75000

winds_GHP_T<-data.frame(rbindlist(mclapply(chrs,function(x){
  tmp<-ghp_t[ghp_t$chrom == x,]
  winds1<-seq(0,max(tmp$pos),by=wind_size)
  winds2<-winds1+wind_size
  
  # Summarise for each
  sum_AF<-data.frame(rbindlist(mclapply(1:length(winds2),function(y){
    tmp2<-tmp[which(tmp$pos<=winds2[y] & tmp$pos>=winds1[y]),]
    if(nrow(tmp2) == 0){
      out <- data.frame("NA")
      out$crit<-sum(tmp2$crit)
    } else {
      out<-data.frame(as.numeric(mean(tmp2$xpehh,na.rm = T)))
      out$crit<-sum(tmp2$crit)
    }
    
    # Tidy
    #out$river<-rownames(out)
    out$normexphh<-as.numeric(mean(tmp2$normxpehh,na.rm = T))
    out$chrom<-x
    out$window<-y
    out$BP1<-as.integer(winds1[y])+1
    out$BP2<-as.integer(winds2[y])
    out$comp<-'IT_GHP'
    colnames(out)<-c('xpehh','crit',"normxpehh",'chrom','window','BP1','BP2','comp')
    out<-out[,c(4:7,1,3,2,8)]
    
    return(out)
  },mc.cores = 3)))
  return(sum_AF)
},mc.cores = 3)))
winds_GHP_T$xpehh<-as.numeric(as.character(winds_GHP_T$xpehh))
winds_GHP_T$ID<-paste0(winds_GHP_T$chrom,":",winds_GHP_T$BP1,"-",winds_GHP_T$BP2)
write.table(winds_GHP_T, "output/selscan/xpehh/GHP_T_norm.xpehh_windows.txt", quote = F, sep = '\t', row.names = F)

winds_GHP_C<-data.frame(rbindlist(mclapply(chrs,function(x){
  tmp<-ghp_c[ghp_c$chrom == x,]
  winds1<-seq(0,max(tmp$pos),by=wind_size)
  winds2<-winds1+wind_size
  
  # Summarise for each
  sum_AF<-data.frame(rbindlist(mclapply(1:length(winds2),function(y){
    tmp2<-tmp[which(tmp$pos<=winds2[y] & tmp$pos>=winds1[y]),]
    if(nrow(tmp2) == 0){
      out <- data.frame("NA")
      out$crit<-sum(tmp2$crit)
    } else {
      out<-data.frame(as.numeric(mean(tmp2$xpehh,na.rm = T)))
      out$crit<-sum(tmp2$crit)
    }
    
    # Tidy
    #out$river<-rownames(out)
    out$normexphh<-as.numeric(mean(tmp2$normxpehh,na.rm = T))
    out$chrom<-x
    out$window<-y
    out$BP1<-as.integer(winds1[y])+1
    out$BP2<-as.integer(winds2[y])
    out$comp<-'IC_GHP'
    colnames(out)<-c('xpehh','crit',"normxpehh",'chrom','window','BP1','BP2','comp')
    out<-out[,c(4:7,1,3,2,8)]
    
    return(out)
  },mc.cores = 3)))
  return(sum_AF)
},mc.cores = 3)))
winds_GHP_C$xpehh<-as.numeric(as.character(winds_GHP_C$xpehh))
winds_GHP_C$ID<-paste0(winds_GHP_C$chrom,":",winds_GHP_C$BP1,"-",winds_GHP_C$BP2)
write.table(winds_GHP_C, "output/selscan/xpehh/GHP_C_norm.xpehh_windows.txt", quote = F, sep = '\t', row.names = F)

winds_GHP_LL<-data.frame(rbindlist(mclapply(chrs,function(x){
  tmp<-ghp_ll[ghp_ll$chrom == x,]
  tmp$pos<-as.numeric(tmp$pos)
  winds1<-seq(0,max(tmp$pos),by=wind_size)
  winds2<-winds1+wind_size
  
  # Summarise for each
  sum_AF<-data.frame(rbindlist(mclapply(1:length(winds2),function(y){
    tmp2<-tmp[which(tmp$pos<=winds2[y] & tmp$pos>=winds1[y]),]
    if(nrow(tmp2) == 0){
      out <- data.frame("NA")
      out$crit<-sum(tmp2$crit)
    } else {
      out<-data.frame(as.numeric(mean(tmp2$xpehh,na.rm = T)))
      out$crit<-sum(tmp2$crit)
    }
    
    # Tidy
    #out$river<-rownames(out)
    out$normexphh<-as.numeric(mean(tmp2$normxpehh,na.rm = T))
    out$chrom<-x
    out$window<-y
    out$BP1<-as.integer(winds1[y])+1
    out$BP2<-as.integer(winds2[y])
    out$comp<-'ILL_GHP'
    colnames(out)<-c('xpehh','crit',"normxpehh",'chrom','window','BP1','BP2','comp')
    out<-out[,c(4:7,1,3,2,8)]
    
    return(out)
  },mc.cores = 3)))
  return(sum_AF)
},mc.cores = 3)))
winds_GHP_LL$xpehh<-as.numeric(as.character(winds_GHP_LL$xpehh))
winds_GHP_LL$ID<-paste0(winds_GHP_LL$chrom,":",winds_GHP_LL$BP1,"-",winds_GHP_LL$BP2)
write.table(winds_GHP_LL, "output/selscan/xpehh/GHP_LL_norm.xpehh_windows.txt", quote = F, sep = '\t', row.names = F)

winds_GHP_UL<-data.frame(rbindlist(mclapply(chrs,function(x){
  tmp<-ghp_ul[ghp_ul$chrom == x,]
  winds1<-seq(0,max(tmp$pos),by=wind_size)
  winds2<-winds1+wind_size
  
  # Summarise for each
  sum_AF<-data.frame(rbindlist(mclapply(1:length(winds2),function(y){
    tmp2<-tmp[which(tmp$pos<=winds2[y] & tmp$pos>=winds1[y]),]
    if(nrow(tmp2) == 0){
      out <- data.frame("NA")
      out$crit<-sum(tmp2$crit)
    } else {
      out<-data.frame(as.numeric(mean(tmp2$xpehh,na.rm = T)))
      out$crit<-sum(tmp2$crit)
    }
    
    # Tidy
    #out$river<-rownames(out)
    out$normexphh<-as.numeric(mean(tmp2$normxpehh,na.rm = T))
    out$chrom<-x
    out$window<-y
    out$BP1<-as.integer(winds1[y])+1
    out$BP2<-as.integer(winds2[y])
    out$comp<-'IUL_GHP'
    colnames(out)<-c('xpehh','crit',"normxpehh",'chrom','window','BP1','BP2','comp')
    out<-out[,c(4:7,1,3,2,8)]
    
    return(out)
  },mc.cores = 3)))
  return(sum_AF)
},mc.cores = 3)))
winds_GHP_UL$xpehh<-as.numeric(as.character(winds_GHP_UL$xpehh))
winds_GHP_UL$ID<-paste0(winds_GHP_UL$chrom,":",winds_GHP_UL$BP1,"-",winds_GHP_UL$BP2)
write.table(winds_GHP_UL, "output/selscan/xpehh/GHP_UL_norm.xpehh_windows.txt", quote = F, sep = '\t', row.names = F)


xp_all<-rbind(winds_GHP_GL,winds_GHP_T,winds_GHP_C,winds_GHP_LL,winds_GHP_UL)
write.table(xp_all, "output/selscan/xpehh/FIBR_xpehh_all.txt", quote = F, sep = '\t', row.names = F)

### Find the outlier windows
# read in the windowed files
c_xp_all<-read.table('output/selscan/xpehh/GHP_C_norm.xpehh_windows.txt',header = T)
t_xp_all<-read.table('output/selscan/xpehh/GHP_T_norm.xpehh_windows.txt',header = T)
ll_xp_all<-read.table('output/selscan/xpehh/GHP_LL_norm.xpehh_windows.txt',header = T)
ul_xp_all<-read.table('output/selscan/xpehh/GHP_UL_norm.xpehh_windows.txt',header = T)

#gl_xp_all[!is.finite(gl_xp_all$normxpehh),]<-NA
c_xp_all[!is.finite(c_xp_all$normxpehh),]<-NA
t_xp_all[!is.finite(t_xp_all$normxpehh),]<-NA
ll_xp_all[!is.finite(ll_xp_all$normxpehh),]<-NA
ul_xp_all[!is.finite(ul_xp_all$normxpehh),]<-NA

#g_out2<-na.omit(gl_xp_all[gl_xp_all$normxpehh>2.5,])
c_out2<-na.omit(c_xp_all[c_xp_all$normxpehh>2.5,])
t_out2<-na.omit(t_xp_all[t_xp_all$normxpehh>2.5,])
ll_out2<-na.omit(ll_xp_all[ll_xp_all$normxpehh>2.5,])
ul_out2<-na.omit(ul_xp_all[ul_xp_all$normxpehh>2.5,])

#write.table(g_out2, "output/selscan/xpehh/outliers/GL_outliers_XP.txt", quote = F, sep = '\t', row.names = F)
write.table(c_out2, "output/selscan/xpehh/outliers/C_outliers_XP.txt", quote = F, sep = '\t', row.names = F)
write.table(t_out2, "output/selscan/xpehh/outliers/T_outliers_XP.txt", quote = F, sep = '\t', row.names = F)
write.table(ll_out2, "output/selscan/xpehh/outliers/LL_outliers_XP.txt", quote = F, sep = '\t', row.names = F)
write.table(ul_out2, "output/selscan/xpehh/outliers/UL_outliers_XP.txt", quote = F, sep = '\t', row.names = F)


