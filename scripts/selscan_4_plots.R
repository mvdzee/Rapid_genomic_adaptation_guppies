library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(data.table)
library(cowplot)
# this script wil produce the plots presented in figure 2 of the paper
# read in the ihh12 file
ih<-read.table("output/selscan/ihh12/FIBR_ihh12_GH_LP.txt",header = T)
ih_chr<-ih[ih$chrom%like%'chr',]
ih_chr<-ih_chr[order(factor(ih_chr$comp,levels=c('ILL_GHP','IUL_GHP','IC_GHP','IT_GHP')),
                     factor(ih_chr$chrom, levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                                                     "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
                                                     "chr20","chr21","chr22","chr23"))),]
pops_ih<-unique(ih_chr$comp)
ih_chr2<-data.frame(matrix(ncol=ncol(ih_chr)+1))
colnames(ih_chr2)<-c(colnames(ih_chr),'gen_wnd')
for(pop in pops_ih){
  tmp<-ih_chr[ih_chr$comp == pop,]
  tmp$gen_wnd<-seq(1,nrow(tmp),1)
  ih_chr2<-rbind(ih_chr2,tmp)
}
ih_chr2<-na.omit(ih_chr2[-1,])
ih_chr2$comp<-factor(ih_chr2$comp,levels=c('ILL_GHP','IUL_GHP','IC_GHP','IT_GHP'))
ih_chr2$chrom<-factor(ih_chr2$chrom, levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                                                "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
                                                "chr20","chr21","chr22","chr23"))

axisdf_ih <- ih_chr2 %>% group_by(chrom) %>% summarize(center=(max(gen_wnd) + min(gen_wnd))/2)
axisdf_ih$chrom<-factor(axisdf_ih$chrom, levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                                                    "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
                                                    "chr20","chr21","chr22","chr23"))
outs_ih<-data.frame(matrix(ncol=ncol(ih_chr2)))
colnames(outs_ih)<-colnames(ih_chr2)
pops_ih<-unique(ih_chr2$comp)
for(pop in pops_ih){
  tmp<-ih_chr2[ih_chr2$comp == pop & ih_chr2$LP_ihh > 5 & ih_chr2$GH_ihh < 5,]
  outs_ih<-rbind(outs_ih,tmp)
}
outs_ih<-outs_ih[-1,]
outs_ih$comp<-factor(outs_ih$comp,levels=c('ILL_GHP','IUL_GHP','IC_GHP','IT_GHP','GLP_GHP'))
outs_ih$chrom<-factor(outs_ih$chrom, levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                                                "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
                                                "chr20","chr21","chr22","chr23"))
outs_ih<-outs_ih[!(outs_ih$comp=="GLP_GHP"),]
col<-c('#686566','grey')

ih_chr2<-ih_chr2[!(ih_chr2$comp=="GLP_GHP"),]
ih_chr2$chrom<-factor(ih_chr2$chrom, levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                                                "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
                                                "chr20","chr21","chr22","chr23"))
ihh_3pop<-c("chr15:5175001-5250000", "chr15:5250001-5325000", "chr15:5325001-5400000", "chr15:5400001-5475000", "chr15:5475001-5550000", "000083F_0.2:150001-225000", "chr9:28500001-28575000")
ihh_2pop<-c( "chr10:25950001-26025000","000104F_0:1200001-1275000","chr1:13800001-13875000","chr9:28575001-28650000,chr1:28050001-28125000","chr13:31575001-31650000","chr15:5100001-5175000","chr4:3525001-3600000,000083F_0.2:150001-225000","chr1:14325001-14400000","chr12:5025001-5100000","chr12:5250001-5325000","chr5:26025001-26100000")
test3<-outs_ih[outs_ih$ID%in%ihh_3pop,]
test2<-outs_ih[outs_ih$ID%in%ihh_2pop,]
outs_ih2<-outs_ih[!(outs_ih$ID%in%test3$ID),]
outs_ih2<-outs_ih2[!(outs_ih2$ID%in%test2$ID),]
outs_ih2$chrom<-factor(outs_ih2$chrom, levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                                                  "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
                                                  "chr20","chr21","chr22","chr23"))

test2$chrom<-factor(test2$chrom, levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                                            "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
                                            "chr20","chr21","chr22","chr23"))
test3$chrom<-factor(test3$chrom, levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                                            "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
                                            "chr20","chr21","chr22","chr23"))
ihh_plot<-ggplot(ih_chr2, aes(x=gen_wnd, y=LP_ihh)) +
  facet_wrap(~comp,nrow=5, strip.position = 'right')+
  geom_point(aes(color=chrom), alpha=0.8, size=1) +
  geom_point(data=outs_ih2,aes(x=gen_wnd,y=LP_ihh),col="#ffdb19",size=1.5,alpha=0.9)+
  geom_point(data=test2,aes(x=gen_wnd,y=LP_ihh),col="#d62976",size=1.5)+
  geom_point(data=test3,aes(x=gen_wnd,y=LP_ihh),col="#fa7e1e",size=1.5)+
  scale_color_manual(values = rep(col, 23)) +
  scale_x_continuous(name="Chromosome", label = seq(1,23,1), breaks= axisdf_ih$center,expand = c(0,0)) +
  scale_y_continuous(name='iHH12',expand = c(0, 0),limits = c(min(ih_chr2$LP_ihh,na.rm = T)-0.5,max(ih_chr2$LP_ihh,na.rm = T)+0.5)) + 
  theme_bw(base_size = 22) +
  theme( 
    plot.title = element_text(hjust = 0.5),
    legend.position="none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size=12),
    strip.text = element_text(size=12),
    axis.title = element_text(size=14)
  )
ihh_plot

xp<-read.table("output/xpehh/FIBR_xpehh_all.txt",header = T)
xp_chr<-xp[xp$chrom%like%'chr',]
xp_chr<-xp_chr[order(xp_chr$comp,factor(xp_chr$chrom, levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                                                                 "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
                                                                 "chr20","chr21","chr22","chr23"))),]
pops<-unique(xp_chr$comp)
xp_chr2<-data.frame(matrix(ncol=ncol(xp_chr)+1))
colnames(xp_chr2)<-c(colnames(xp_chr),'gen_wnd')
for(pop in pops){
  tmp<-xp_chr[xp_chr$comp == pop,]
  tmp$gen_wnd<-seq(1,nrow(tmp),1)
  xp_chr2<-rbind(xp_chr2,tmp)
}
xp_chr2<-na.omit(xp_chr2[-1,])
xp_chr2$comp<-factor(xp_chr2$comp,levels=c('ILL_GHP','IUL_GHP','IC_GHP','IT_GHP','GLP_GHP'))
xp_chr2$chrom<-factor(xp_chr2$chrom, levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                                                "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
                                                "chr20","chr21","chr22","chr23"))

axisdf <- xp_chr2 %>% group_by(chrom) %>% summarize(center=( max(gen_wnd) + min(gen_wnd) ) / 2 )
axisdf$chrom<-factor(axisdf$chrom, levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                                              "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
                                              "chr20","chr21","chr22","chr23"))
outs<-data.frame(matrix(ncol=ncol(xp_chr2)))
colnames(outs)<-colnames(xp_chr2)
pops<-unique(xp_chr2$comp)
for(pop in pops){
  tmp<-xp_chr2[xp_chr2$comp == pop & xp_chr2$normxpehh > 2.5,]
  outs<-rbind(outs,tmp)
}
outs<-outs[-1,]
outs$comp<-factor(outs$comp,levels=c('ILL_GHP','IUL_GHP','IC_GHP','IT_GHP','GLP_GHP'))
outs$chrom<-factor(outs$chrom, levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                                          "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
                                          "chr20","chr21","chr22","chr23"))
outs<-outs[!(outs$comp=="GLP_GHP"),]
xp_chr2<-xp_chr2[!(xp_chr2$comp=="GLP_GHP"),]

col<-c('#686566','grey')

xp_2pop<-c("chr15:5250001-5325000","chr15:5325001-5400000","chr15:5400001-5475000","chr15:5475001-5550000","chr18:21900001-21975000","chr15:5100001-5175000","chr15:5550001-5625000","chr4:3600001-3675000","chr4:3675001-3750000","chr15:5025001-5100000","chr5:26175001-26250000","chr5:26250001-26325000")
test2xp<-outs[outs$ID%in%xp_2pop,]
outs2<-outs[!(outs$ID%in%test2xp$ID),]
outs2$chrom<-factor(outs2$chrom, levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                                            "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19",
                                            "chr20","chr21","chr22","chr23"))
xp_plot<-ggplot(xp_chr2, aes(x=gen_wnd, y=normxpehh)) +
  # Show all points
  facet_wrap(~comp,nrow=5, strip.position = 'right')+
  geom_point(aes(color=chrom), alpha=0.8, size=1) +
  geom_point(data=outs2,aes(x=gen_wnd,y=normxpehh),col="#ffdb19",size=1.5,alpha=0.9)+
  geom_point(data=test2xp,aes(x=gen_wnd,y=normxpehh),col="#d62976",size=1.5)+
  scale_color_manual(values = rep(col, 23)) +
  scale_x_continuous(name="",label = seq(1,23,1), breaks= axisdf$center,expand = c(0,0)) +
  scale_y_continuous(name='XP-EHH',expand = c(0, 0),limits = c(min(xp_chr2$normxpehh,na.rm = T)-0.5,max(xp_chr2$normxpehh,na.rm = T)+0.5)) +
  geom_hline(yintercept = 0,col='white') +
  theme_bw(base_size = 22) +
  theme( 
    plot.title = element_text(hjust = 0.5),
    legend.position="none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size=12),
    strip.text = element_text(size=12),
    axis.title = element_text(size=14)
  )
xp_plot

pdf('figures/selscan/xp_ih_whole_genome_xp2.5_ihh5.pdf',width = 8, height=11)
plot_grid(xp_plot,ihh_plot,ncol=1, align = 'v',labels='AUTO',label_size = 14)
dev.off()

