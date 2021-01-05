library(ggplot2)
library(grid)
library(plyr)
library(ggpubr)
library(gridExtra)
library(dplyr)
library(data.table)
library(viridis)
library(cowplot)

### This script will process the ROH and create the plot seen in fig 1C
### Read in the data
c_hom<-read.table("data/C_500K_1het.hom", header=T)
t_hom<-read.table("data/T_500K_1het.hom", header=TRUE)
#glp_hom<-read.table("data/ROH/GLP_500K_1het.hom", header=TRUE)
ghp_hom<-read.table("data/GHP_500K_1het.hom", header=TRUE)
ll_hom<-read.table("data/ILL_500K_1het.hom", header=TRUE)
ul_hom<-read.table("data/IUL_500K_1het.hom", header=TRUE)

# remove the sex chromosome
c_hom<-c_hom[!(c_hom$CHR==12),]
t_hom<-t_hom[!(t_hom$CHR==12),]
ll_hom<-ll_hom[!(ll_hom$CHR==12),]
ul_hom<-ul_hom[!(ul_hom$CHR==12),]
#glp_hom<-glp_hom[!(glp_hom$CHR==12),]
ghp_hom<-ghp_hom[!(ghp_hom$CHR==12),]

### Add population group
c_hom$POP<-as.factor(rep('IC',nrow(c_hom)))
t_hom$POP<-as.factor(rep('IT',nrow(t_hom)))
#glp_hom$POP<-as.factor(rep('GLP',nrow(glp_hom)))
ghp_hom$POP<-as.factor(rep('GHP',nrow(ghp_hom)))
ll_hom$POP<-as.factor(rep('ILL',nrow(ll_hom)))
ul_hom$POP<-as.factor(rep('IUL',nrow(ul_hom)))

## Add MB column
c_hom$MB<-c_hom$KB/1000
t_hom$MB<-t_hom$KB/1000
#glp_hom$MB<-glp_hom$KB/1000
ghp_hom$MB<-ghp_hom$KB/1000
ll_hom$MB<-ll_hom$KB/1000
ul_hom$MB<-ul_hom$KB/1000

### Bin the data in 4 bins, and summarise: 1) the total length of ROHs in a bin and 2) the number of ROH's in a bin 
c_hom$GROUP<-cut(c_hom$MB, breaks=c(0.5,0.75,1,1.5,2,2.5,Inf),labels=c("0.50-0.75","0.75-1.0","1.0-1.5","1.5-2.0","2.0-2.5",">2.5"))
# sum per group per individual
c_sumind<-aggregate(list(Mb=c_hom$MB), by=list(group=c_hom$GROUP,ind=c_hom$IID), FUN=sum)
# average length per group
c_sum<-aggregate(list(sum=c_sumind$Mb), by=list(Category=c_sumind$group), FUN=mean)
c_sum$POP<-rep('IC',nrow(c_sum))
c_sum$prop<-c_sum$sum/sum(c_sum$sum)
c_sum$propBUCK<-c_sum$sum/739.3
# counts per group per individual
c_freq2<- ddply(c_hom, .(c_hom$GROUP, c_hom$IID), nrow)
# average count per group
c_freq<- aggregate(list(mean_count=c_freq2$V1), by=list(Category=c_freq2$`c_hom$GROUP`), FUN=mean)
c_freq$POP<-rep('IC',nrow(c_freq))
c_freq$prop<-c_freq$mean_count/sum(c_freq$mean_count)

t_hom$GROUP<-cut(t_hom$MB, breaks=c(0.5,0.75,1,1.5,2,2.5,Inf),labels=c("0.50-0.75","0.75-1.0","1.0-1.5","1.5-2.0","2.0-2.5",">2.5"))
t_sumind<-aggregate(list(MB=t_hom$MB), by=list(group=t_hom$GROUP,ind=t_hom$IID), FUN=sum)
t_sum<-aggregate(list(sum=t_sumind$MB), by=list(Category=t_sumind$group), FUN=mean)
t_sum$POP<-rep('IT',nrow(t_sum))
t_sum$prop<-t_sum$sum/sum(t_sum$sum)
t_sum$propBUCK<-t_sum$sum/739.3
t_freq2<- ddply(t_hom, .(t_hom$GROUP, t_hom$IID), nrow)
t_freq<- aggregate(list(mean_count=t_freq2$V1), by=list(Category=t_freq2$`t_hom$GROUP`), FUN=mean)
t_freq$POP<-rep('IT',nrow(t_freq))
t_freq$prop<-t_freq$mean_count/sum(t_freq$mean_count)

ghp_hom$GROUP<-cut(ghp_hom$MB, breaks=c(0.5,0.75,1,1.5,2,2.5,Inf),labels=c("0.50-0.75","0.75-1.0","1.0-1.5","1.5-2.0","2.0-2.5",">2.5"))
ghp_sumind<-aggregate(list(MB=ghp_hom$MB), by=list(group=ghp_hom$GROUP,ind=ghp_hom$IID), FUN=sum)
ghp_sum<-aggregate(list(sum=ghp_sumind$MB), by=list(Category=ghp_sumind$group), FUN=mean)
ghp_sum$POP<-rep('GHP',nrow(ghp_sum))
ghp_sum$prop<-ghp_sum$sum/sum(ghp_sum$sum)
ghp_sum$propBUCK<-ghp_sum$sum/739.3
ghp_freq2<- ddply(ghp_hom, .(ghp_hom$GROUP, ghp_hom$IID), nrow)
ghp_freq<- aggregate(list(mean_count=ghp_freq2$V1), by=list(Category=ghp_freq2$`ghp_hom$GROUP`), FUN=mean)
ghp_freq$POP<-rep('GHP',nrow(ghp_freq))
ghp_freq$prop<-ghp_freq$mean_count/sum(ghp_freq$mean_count)

ll_hom$GROUP<-cut(ll_hom$MB, breaks=c(0.5,0.75,1,1.5,2,2.5,Inf),labels=c("0.50-0.75","0.75-1.0","1.0-1.5","1.5-2.0","2.0-2.5",">2.5"))
ll_sumind<-aggregate(list(MB=ll_hom$MB), by=list(group=ll_hom$GROUP,ind=ll_hom$IID), FUN=sum)
ll_sum<-aggregate(list(sum=ll_sumind$MB), by=list(Category=ll_sumind$group), FUN=mean)
ll_sum$POP<-rep('ILL',nrow(ll_sum))
ll_sum$prop<-ll_sum$sum/sum(ll_sum$sum)
ll_sum$propBUCK<-ll_sum$sum/739.3
ll_freq2<- ddply(ll_hom, .(ll_hom$GROUP, ll_hom$IID), nrow)
ll_freq<- aggregate(list(mean_count=ll_freq2$V1), by=list(Category=ll_freq2$`ll_hom$GROUP`), FUN=mean)
ll_freq$POP<-rep('ILL',nrow(ll_freq))
ll_freq$prop<-ll_freq$mean_count/sum(ll_freq$mean_count)

ul_hom$GROUP<-cut(ul_hom$MB, breaks=c(0.5,0.75,1,1.5,2,2.5,Inf),labels=c("0.50-0.75","0.75-1.0","1.0-1.5","1.5-2.0","2.0-2.5",">2.5"))
ul_sumind<-aggregate(list(MB=ul_hom$MB), by=list(group=ul_hom$GROUP,ind=ul_hom$IID), FUN=sum)
ul_sum<-aggregate(list(sum=ul_sumind$MB), by=list(Category=ul_sumind$group), FUN=mean)
ul_sum$POP<-rep('IUL',nrow(ul_sum))
ul_sum$prop<-ul_sum$sum/sum(ul_sum$sum)
ul_sum$propBUCK<-ul_sum$sum/739.3
ul_freq2<- ddply(ul_hom, .(ul_hom$GROUP, ul_hom$IID), nrow)
ul_freq<- aggregate(list(mean_count=ul_freq2$V1), by=list(Category=ul_freq2$`ul_hom$GROUP`), FUN=mean)
ul_freq$POP<-rep('IUL',nrow(ul_freq))
ul_freq$prop<-ul_freq$mean_count/sum(ul_freq$mean_count)

freq_all_gl<-rbind(ghp_freq,t_freq,c_freq,ul_freq,ll_freq)
freq_all_gl$POP<-factor(freq_all_gl$POP, levels = c('ILL','IUL','IC','IT','GHP'))
sum_all_gl<-rbind(ghp_sum,t_sum,c_sum,ul_sum,ll_sum)
sum_all_gl$POP<-factor(sum_all_gl$POP, levels = c('ILL','IUL','IC','IT','GHP'))

write.table(freq_all_gl,"output/ROH/freq_of_ROH_auto_500K.txt", quote = F, sep = '\t', row.names = F)
write.table(sum_all_gl,"output/ROH/sum_of_ROH_auto_500K.txt", quote = F, sep = '\t', row.names = F)

### plot bars. On x is the size bins. Y = the sum of the ROH's for that bin in each population
png('figures/ROH/sum_of_ROH_auto_500K_bars.png',width = 800,height = 300)
sums<-ggplot(data=sum_all_gl, aes(x=POP,y=sum,fill=Category))+
  geom_bar(stat="identity",width = 0.9,position=position_dodge2(width = 2, preserve = "single"),color='#238443')+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size=20),
        axis.title=element_text(size=20),
        plot.title = element_text(hjust = 0.5,size=20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size=20),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
        legend.position = 'none') +
  scale_fill_manual(values = c('#f7fcb9','#d9f0a3','#addd8e','#78c679','#41ab5d','#238443','#005a32'))+
  scale_y_continuous(name='Sum of ROHs (Mb)',expand = c(0,0))+
  xlab('Population')+
  labs(fill = "Size class (Mb)")+
  coord_cartesian(ylim=c(0,15))
sums
dev.off()

# plot the nr of ROH
png('figures/ROH/freqs_of_ROH_auto_500K_bars_line.png',width = 800,height = 300)
freqs<-ggplot(data=freq_all_gl, aes(x=POP,y=mean_count,fill=Category))+
  geom_bar(stat="identity",width = 0.9,position=position_dodge2(width = 2, preserve = "single"),color='#238443')+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size=20),
        axis.title = element_text(size=20),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=20),
        plot.title = element_text(hjust = 0.5,size=20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size=20),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
        legend.position = 'none') +
  scale_fill_manual(values = c('#f7fcb9','#d9f0a3','#addd8e','#78c679','#41ab5d','#238443','#005a32'))+
  scale_y_continuous(name='Count of ROHs (Mb)',expand = c(0,0))+
  xlab('Population')+
  labs(fill = "Size class (Mb)")+
  coord_cartesian(ylim=c(0,20))
freqs
dev.off()

# get the legend
legend_b <- get_legend(
  freqs + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom",
          legend.box.margin = margin(0, 10, 0, 20))
)

pdf('figures/FIBR_ROH_final_lines.pdf',width=8,height=8)
plot_grid(freqs,sums,legend_b, rel_heights = c(1,1,0.2),ncol=1)
dev.off()

### table as in supplementary table S4
# total ROH per pop
nrow(c_hom)
nrow(t_hom)
nrow(ll_hom)
nrow(ul_hom)
nrow(ghp_hom)

# mean Mb in ROH, min length + max length
c_length<-aggregate(c_hom$MB, by=list(ind=c_hom$IID), FUN=sum)
c_mean<-mean(c_length$x)
min(c_hom$MB)
max(c_hom$MB)

t_length<-aggregate(t_hom$MB, by=list(ind=t_hom$IID), FUN=sum)
eh_mean<-mean(t_length$x)
min(t_hom$MB)
max(t_hom$MB)

ll_length<-aggregate(ll_hom$MB, by=list(ind=ll_hom$IID), FUN=sum)
tl_mean<-mean(ll_length$x)
min(ll_hom$MB)
max(ll_hom$MB)

ul_length<-aggregate(ul_hom$MB, by=list(ind=ul_hom$IID), FUN=sum)
th_mean<-mean(ul_length$x)
min(ul_hom$MB)
max(ul_hom$MB)

ghp_length<-aggregate(ghp_hom$MB, by=list(ind=ghp_hom$IID), FUN=sum)
gh_mean<-mean(ghp_length$x)
min(ghp_hom$MB)
max(ghp_hom$MB)

# calculate the inbreeding coefficient
## Froh = SUM(LRroh)/Lauto
c_froh<-aggregate(list(Mb=c_hom$MB), by=list(ind=c_hom$IID), FUN=sum)
c_froh$FROH<-c_froh$Mb/739.3
mean(c_froh$FROH)
median(c_froh$FROH)
min(c_froh$FROH)
max(c_froh$FROH)

t_froh<-aggregate(list(Mb=t_hom$MB), by=list(ind=t_hom$IID), FUN=sum)
t_froh$FROH<-t_froh$Mb/739.3
mean(t_froh$FROH)
median(t_froh$FROH)
min(t_froh$FROH)
max(t_froh$FROH)

ghp_froh<-aggregate(list(Mb=ghp_hom$MB), by=list(ind=ghp_hom$IID), FUN=sum)
ghp_froh$FROH<-ghp_froh$Mb/739.3
mean(ghp_froh$FROH)
median(ghp_froh$FROH)
min(ghp_froh$FROH)
max(ghp_froh$FROH)

glp_froh<-aggregate(list(Mb=glp_hom$MB), by=list(ind=glp_hom$IID), FUN=sum)
glp_froh$FROH<-glp_froh$Mb/739.3
median(glp_froh$FROH)
min(glp_froh$FROH)
max(glp_froh$FROH)

ll_froh<-aggregate(list(Mb=ll_hom$MB), by=list(ind=ll_hom$IID), FUN=sum)
ll_froh$FROH<-ll_froh$Mb/739.3
median(ll_froh$FROH)
min(ll_froh$FROH)
max(ll_froh$FROH)

ul_froh<-aggregate(list(Mb=ul_hom$MB), by=list(ind=ul_hom$IID), FUN=sum)
ul_froh$FROH<-ul_froh$Mb/739.3
median(ul_froh$FROH)
min(ul_froh$FROH)
max(ul_froh$FROH)
