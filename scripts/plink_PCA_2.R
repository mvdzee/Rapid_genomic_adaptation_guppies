# this R script creates the PCA plot in figure 1B

library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)

# read in the eigenvectors
eigenvec_table <- read.table('data/PCA/FIBR_PCA_noGL.eigenvec', header = TRUE)
eigenvec_table <- eigenvec_table[-2]
eigenvec_table <- eigenvec_table[-1]
# provide the population names in the order they appear in your vcf (alphabetical)
eigenvec_table$Populations <- as.character(c(rep("IC (intro)", 15), rep("Guanapo (HP)", 19), rep("ILL (intro)", 13), rep("IT (intro)", 14), rep("IUL (intro)", 15)))

# set the colour scheme per population
palette<-(c(rep("#088da5"),rep("#dd3497"),rep("#fa9fb5"),rep("#756bb1"),rep("#bcbddc")))

# read in the eigenvalues and calculate the percentages
eigenval <- read.table('data/PCA/FIBR_PCA_noGL.eigenval', header = F)
percentage <- round(eigenval$V1/sum(eigenval$V1)*100,2)
percentage <- paste(colnames(eigenvec_table),"(",paste(as.character(percentage),"%",")",sep=""))

pop_order <- c("Guanapo (HP)","ILL (intro)","IUL (intro)","IC (intro)","IT (intro)")
raw_levels<- c("Guanapo (HP)","ILL (intro)","IUL (intro)","IC (intro)","IT (intro)")

for (i in 1:length(pop_order)){
  eigenvec_table[eigenvec_table$Populations==raw_levels[i],length(colnames(eigenvec_table))] <- pop_order[i]
}

eigenvec_table$Populations <- factor(eigenvec_table$Populations, levels=pop_order)

# plot PC1 & PC2
pdf('figures/PCA/FIBR_noGL_pc1_pc2_legend.pdf',width=20, height=8)
p1<-ggplot(eigenvec_table,aes(x=PC1,y=PC2,color=Populations)) +
  geom_point(size = 4) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.title=element_text(size=20),
        legend.text = element_text(size=20),
        axis.text = element_text(size=20),
        axis.title = element_text(size=20)) +
  scale_colour_manual(values=c('Guanapo (HP)'="#088da5",'IC (intro)'="#756BB1",'IT (intro)'="#bcbddc",'ILL (intro)'="#dd3497",'IUL (intro)'="#fa9fb5"))+
  stat_ellipse(level=0.95,show.legend = F) +
  xlab(percentage[1]) +
  ylab(percentage[2]) +
  guides(fill=guide_legend(title="Populations",size=40), color=guide_legend(override.aes=list(size=8)))
print(p1)
dev.off()

# plot PC1 & PC3
pdf('figures/PCA/FIBR_noGL_pc1_pc3.pdf',width=8, height=6)
p2<-ggplot(eigenvec_table,aes(x=PC1,y=PC3,color=Populations)) +
  geom_point(size = 4) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.title=element_text(size=20),
        legend.text = element_text(size=20),
        axis.text = element_text(size=20),
        axis.title = element_text(size=20)) +
  scale_colour_manual(values=c('Guanapo (HP)'="#088da5",'IC (intro)'="#756BB1",'IT (intro)'="#bcbddc",'ILL (intro)'="#dd3497",'IUL (intro)'="#fa9fb5"))+
  stat_ellipse(level=0.95,show.legend = F) +
  xlab(percentage[1]) +
  ylab(percentage[3]) +
  guides(fill=guide_legend(title="Populations"), color=guide_legend(override.aes=list(size=4)))
print(p2)
dev.off()

# plot PC2 & PC3
pdf('figures/PCA/FIBR_noGL_pc2_pc3.pdf',width=8, height=6)
p3<-ggplot(eigenvec_table,aes(x=PC2,y=PC3,color=Populations)) +
  geom_point(size = 4) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.title=element_text(size=20),
        legend.text = element_text(size=20),
        axis.text = element_text(size=20),
        axis.title = element_text(size=20)) +
  scale_colour_manual(values=c('Guanapo (HP)'="#088da5",'IC (intro)'="#756BB1",'IT (intro)'="#bcbddc",'ILL (intro)'="#dd3497",'IUL (intro)'="#fa9fb5"))+
  stat_ellipse(level=0.95,show.legend = F) +
  xlab(percentage[2]) +
  ylab(percentage[3]) +
  guides(fill=guide_legend(title="Populations"), color=guide_legend(override.aes=list(size=4)))
print(p3)
dev.off()

# plot together which ever ones you need
prow <- plot_grid(
  p1 + theme(legend.position="none"),
  #p2 + theme(legend.position="none"),
  #p3 + theme(legend.position="none"),
  align = 'vh',
  #labels = c("A", "B", "C"),
  label_size = 20,
  hjust = -1,
  nrow = 1
)
prow

legend_b <- get_legend(
  p1 +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)

# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
pdf('figures/PCA/FIBR_noGL_PC12.pdf',width=12, height=4.5,onefile=FALSE)
plot_grid(prow, legend_b,ncol=1, rel_heights = c(3, .4))
dev.off()
