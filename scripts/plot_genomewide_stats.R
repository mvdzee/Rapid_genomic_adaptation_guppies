# this R script will plot the distributions of Expected heterozygosity, nucleotide diversity and Tajima's D for all populations
# as in fig 1D
library(ggplot2)
library(ggpubr)
library(cowplot)

# read in the files
het<-na.omit(read.table("output/heterozygosity/heterozygosity_FIBR.txt",header = T))
pi<-read.table("output/popgenome/FIBR_Pi_processed.txt",header = T)
td<-read.table("output/popgenome/FIBR_TD_processed.txt",header = T)

# remove any populations you don't want to include
het2<-het[!(het$comp=='GLP'),]

# plot heterozygosity
het_plot<-ggplot(data=het2,aes(x=exp_het)) +
  geom_freqpoly(aes(color=comp),binwidth = 0.009,size=1.5,alpha=0.8)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=20),
        axis.title = element_text(size=20),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(size=16),
        legend.text = element_text(size=16),
        plot.margin = margin(c(0.5,0.7,0.1,0.1),unit = "cm"))+
  scale_colour_manual(values=c(GHP="#088da5",GLP="#9cd1db",IC="#756BB1",IT="#bcbddc",ILL="#dd3497",IUL="#fa9fb5"))+
  xlab(expression(H[e]))+
  scale_x_continuous(expand = c(0, 0),limits = c(0,0.6),breaks = seq(0,0.6,0.2),labels=c("0","0.2","0.4","0.6"))+
  scale_y_continuous(name="Count",expand = c(0, 0),limits=c(0,700),breaks=seq(0,800,200))
het_plot

# plot nucleotide diversity
pi<-pi[!(pi$comp=="GLP"),]
pi_plot<-ggplot(data=pi,aes(x=Pi)) +
  geom_freqpoly(aes(color=comp),binwidth = 0.0001,size=1.5,alpha=0.8)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=20),
        axis.title = element_text(size=20),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "None",
        plot.margin = margin(c(0.5,0.7,0.1,0.1),unit = "cm"))+
  scale_colour_manual(values=c(GHP="#088da5",GLP="#9cd1db",IC="#756BB1",IT="#bcbddc",ILL="#dd3497",IUL="#fa9fb5"))+
  xlab(expression(pi))+
  scale_x_continuous(expand = c(0, 0),limits=c(0,0.01),breaks = seq(0,0.01,0.005),labels = seq(0,0.01,0.005))+
  scale_y_continuous(expand = c(0, 0),limits=c(0,600),breaks=seq(0,600,200))
pi_plot

# plot tajima's D
td_g<-td[!(td$comp=="GLP"),]
td_plot<-ggplot(data=td_g,aes(x=TD)) +
  geom_freqpoly(aes(color=comp),binwidth = 0.1,size=1.5,alpha=0.8)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=20),
        axis.title = element_text(size=20),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "None",
        plot.margin = margin(c(0.5,0.7,0.1,0.1),unit = "cm"))+
  scale_colour_manual(values=c(GHP="#088da5",GLP="#9cd1db",IC="#756BB1",IT="#bcbddc",ILL="#dd3497",IUL="#fa9fb5"))+
  xlab("Tajima's D")+
  scale_x_continuous(expand = c(0, 0),limits = c(-3,5),breaks = seq(-4,6,2)) +
  scale_y_continuous(name="Count",expand = c(0, 0),limits=c(0,600),breaks=seq(0,800,100))
td_plot

# combine the 3 plots
prow <- plot_grid(
  het_plot + theme(legend.position="None"),
  pi_plot + theme(legend.position="None"),
  td_plot + theme(legend.position="None"),
  align = 'vh',
  #labels = c("a", "b", "c"),
  #label_size = 18,
  nrow = 1
)
prow

legend_b <- get_legend(
  het_plot +
    guides(color = guide_legend(nrow = 1,override.aes = list(size = 2),title="Population")) +
    theme(legend.position = "bottom",legend.text=element_text(size=14), legend.title = element_text(size=16))
)
  
pdf('figures/heterozygosity/FIBR_general_stats_panel_nolegend.pdf',width=13.5, height=4.5)
plot_grid(prow,ncol=1, rel_heights = c(3, .4))
dev.off()


