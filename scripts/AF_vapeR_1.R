# Genome-wide Test of Eigen analysis
source("scripts/AF_vapeR_functions.R")
library(tidyverse)
library(ggplot2)
# Test inputs
test_vcf <- "data/FIBR_finalVCF.vcf.gz"
window_snps <- 200
CORES <- 2
permutations <- 1000
# Make the vectors
test_vectors=list(c("GH","LL"),c("GH","UL"),c("GH","C"),c("GH","T"))
names(test_vectors) <- c("ILL","IUL","IC","IT")
# Run over the following chr
chrs <- c(paste0("chr",1:23))
chr_res <- lapply(chrs,function(chr){
  
  message(paste0("STARTING ",chr))
  
  # Make the input
  system(paste0("bcftools view -r ",chr," ",test_vcf," > output/tmp.vcf"))
  vcf_in <- read.vcfR("output/AF_vapeR/tmp.vcf")
  system("rm -f output/AF_vapeR/tmp.vcf")
  
  # Make the popmaps
  test_popmap <- data.frame(inds=colnames(vcf_in@gt)[2:ncol(vcf_in@gt)])
  test_popmap$pop <- gsub('[0-9]+', '', test_popmap$inds)
  test_popmap$pop <- gsub('_F', '', test_popmap$pop)
  test_popmap$pop <- gsub('_M', '', test_popmap$pop)
  
  # Calculate allele frequency matrices
  AF_input <- calc_AF_vectors(vcf = vcf_in,
                              window_size = window_snps,
                              popmap = test_popmap,
                              vectors = test_vectors,
                              n_cores = CORES)
  
  # # Calculate null frequencies
  null_input <- calc_null_AF_vectors(vcf = vcf_in,
                                     window_size = window_snps,
                                     popmap = test_popmap,
                                     vectors = test_vectors,
                                     n_cores = CORES,
                                     perm = permutations)
  
  # Calculate Eigen statistics 
  eigen_res <- lapply(AF_input,eigen_analyse_vectors)
  
  # Get null matrix
  null_cutoffs <- find_null_cutoff(null_res = null_input,
                                   cutoffs = c(0.95,0.99,0.999))
  
  # Save both
  saveRDS(list(AF_input,null_input,eigen_res,null_cutoffs),
          paste0("output/AF_vapeR/FIBR_",chr,"_AF_eigen_res_windsize_",window_snps,".rds"))
  
  # Return
  # return(list(AF_input,eigen_res,null_cutoffs))
  return(list(AF_input,eigen_res))
})

# Read in the rds results if needs be...
chr_res <- lapply(chrs,function(chr){
  return(readRDS(paste0("output/AF_vapeR/FIBR_",chr,"_AF_eigen_res_windsize_",window_snps,".rds")))
})
# Just get the chr AF_inputs
AF_input_chrs <-  lapply(chr_res,function(x){return(x[[1]])})
names(AF_input_chrs) <- chrs
# Pull all the eigen_res together
eigen_res_list_test <- lapply(chr_res,function(x){return(x[[3]])})
names(eigen_res_list_test) <- chrs
# Collate the null cutoffs
null_cutoff_list <- lapply(chr_res,function(x){return(x[[4]][,1])})
names(null_cutoff_list) <- chrs
#null_cutoff_list$`000094F_0` <- null_cutoff_list$chr20
# Make a single null cutoff
null_input_list <- lapply(chr_res,function(x){return(x[[2]])})
all_nulls <- unlist(null_input_list, recursive=FALSE)
all_null_cutoffs <- find_null_cutoff(null_res = all_nulls,
                                     cutoffs = c(0.95,0.99,0.999,0.9999))


# Visualise the null cutoffs...
null_eigens <- lapply(all_nulls,eigen_analyse_vectors)
sum_nulls <- sum_eigenvals(null_eigens)
colnames(sum_nulls) <- paste0("Eigenvector ",1:ncol(sum_nulls))
reshape2::melt(sum_nulls) %>%
  ggplot(aes(x=value,fill=Var2))+
  geom_density(alpha=0.5)+
  facet_wrap(~Var2,scales="free_y",ncol=1)+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text.y=element_blank(),
        axis.title = element_text(size=16),
        axis.text.x = element_text(size=15),
        strip.text = element_text(size=14))+
  labs(x="Eigenvalue",y="Density")

# Plot Genome-wide
test_list <- null_cutoff_list
genome_figs <- eigenval_plot_genome(eigen_res_list=eigen_res_list_test,
                                    cutoffs=all_null_cutoffs[,3])
#genome_figs$`Eigenvector 1`
#genome_figs$`Eigenvector 2`
#genome_figs$`Eigenvector 3`
#genome_figs$`Eigenvector 4`
# Plot together
pdf("figures/AF_vapeR/FIBR_genome_99.9CI_200SNPs.pdf",width=12,height=6)
cowplot::plot_grid(genome_figs[[1]] +theme(strip.text = element_blank(),
                                           axis.title.x = element_blank(),plot.title = element_text(size=16),axis.ticks.y=element_line()),
                   genome_figs[[2]] + theme(plot.title = element_text(size=16),axis.ticks.y=element_line()) ,nrow = 2,align = 'v',rel_heights = c(5,6))
dev.off()

# Plot per chr...
chr_plots <- lapply(chrs,function(chr){
  fig <- eigenval_plot(eigen_res = eigen_res_list_test[[chr]],
                       cutoffs=all_null_cutoffs[,3])
  return(fig)
})

# Fetch all our significant windows...
chr_signif <- lapply(chrs,function(chr){
  # Fetch...
  signif <- signif_eigen_windows(eigen_res_list_test[[chr]],cutoffs = all_null_cutoffs[,3])
  return(signif)
})

# We want to look at eigenvector 2
eig1_signif <- na.omit(unlist(lapply(chr_signif,function(chr){return(chr[[1]])})))
eig2_signif <- na.omit(unlist(lapply(chr_signif,function(chr){return(chr[[2]])})))
write.table(eig1_signif, "output/AF_vapeR/FIBR_eig1_outliers_99.9_200SNPs.txt", quote = F, sep = '\t', row.names = F)

# Summarise eigenvector 1 windows...
FIBR_eigenvec1_window_summaries <- data.frame(rbindlist(lapply(eig1_signif,function(wind){
  # First find the necessary chromosome
  eigen_chr <- strsplit(wind,":")[[1]][1]
  # Merge all the eigen_res
  summarise_window_parallelism(wind,eigen_res_list_test[[eigen_chr]],loading_cutoff = 0.3,eigenvector = 1)
})))

# Summarise eigenvector 2 windows...
FIBR_eigenvec2_window_summaries <- data.frame(rbindlist(lapply(eig2_signif,function(wind){
  # First find the necessary chromosome
  eigen_chr <- strsplit(wind,":")[[1]][1]
  # Merge all the eigen_res
  tmp1 <- summarise_window_parallelism(wind,eigen_res_list_test[[eigen_chr]],loading_cutoff = 0.3,eigenvector = 1)
  tmp2 <- summarise_window_parallelism(wind,eigen_res_list_test[[eigen_chr]],loading_cutoff = 0.3,eigenvector = 2)
  return(rbind(tmp1,tmp2))
})))

# Examine A matrices...
eig1_A <- lapply(eig1_signif,function(x){
  # First find the necessary chromosome
  eigen_chr <- strsplit(x,":")[[1]][1]
  # Fetch A matrix
  mat <- data.frame(eigen_res_list_test[[eigen_chr]][[x]][[3]])
  colnames(mat) <- paste0("Eig",1:ncol(mat))
  mat$pos <- as.integer(rownames(mat))
  return(mat)
})
names(eig1_A) <- eig1_signif
eig2_A <- lapply(eig2_signif,function(x){
  # First find the necessary chromosome
  eigen_chr <- strsplit(x,":")[[1]][1]
  # Fetch A matrix
  mat <- data.frame(eigen_res_list_test[[eigen_chr]][[x]][[3]])
  colnames(mat) <- paste0("Eig",1:ncol(mat))
  mat$pos <- as.integer(rownames(mat))
  return(eigen_res_list_test[[eigen_chr]][[x]][[3]])
})
names(eig2_A) <- eig2_signif
#################################################################





