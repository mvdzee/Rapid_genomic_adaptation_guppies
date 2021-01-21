# Dependencies
library(vcfR)
library(data.table)
library(parallel)
library(pegas)
library(ggplot2)
library(adegenet)

################################################################################
calc_AF_vectors <- function(vcf=NULL, # Object of type vcfR, one chrom
                            window_size=100, # window size in SNPs
                            popmap=NULL, # Popmap
                            vectors=NULL, # List of all vectors from ancestor to descendent
                            n_cores=1,
                            normalise=TRUE,
                            end_cutoff=20){
  
  # Check for only one chr and make SNP windows
  if(length(unique(vcf@fix[,1])) > 1){
    stop("VCF contains more than one chromosome...")
  }
  winds <- seq(1,nrow(vcf@fix),window_size)
  winds2 <- winds+window_size-1
  winds2[length(winds2)] <- nrow(vcf@fix)
  
  # Trim any windows that are less than a user-defined cutoff SNPs...
  if(max(winds2)-max(winds) < end_cutoff){
    winds <- winds[1:(length(winds)-1)]
    winds2 <- winds2[1:(length(winds2)-1)]
  }
  
  # How many windows to estimate
  message(paste0("Calculating AF vectors for ",length(winds)," windows with ",window_size," SNPs each"))
  
  # Run in parallel over all windows
  window_list <- mclapply(1:length(winds),function(x){
    
    # Make a "row" in the matrix for each vector
    vector_out <- lapply(vectors,function(vec){
      
      # Subset the VCF for each pop in the vector...
      tmp_vcf <- vcf[winds[x]:winds2[x],c(1,
                                          which(colnames(vcf@gt) %in% popmap[popmap[,2] == vec[1],1]),
                                          which(colnames(vcf@gt) %in% popmap[popmap[,2] == vec[2],1]))]
      
      # Calculate the allele frequencies and normalise
      af <- by(vcfR2loci(tmp_vcf),c(popmap[popmap[,2] == vec[1],2],
                                    popmap[popmap[,2] == vec[2],2]))
      
      # Get AF change
      af_vals <- sapply(af,function(i){
        # Mask bum windows
        if(ncol(i) == 0){
          return(0)
        } else {
          rownames(i) <- sort(vec)
          freq1 <- i[vec[1],1]/sum(i[vec[1],])
          freq2 <- i[vec[2],1]/sum(i[vec[2],])
          return(freq2 - freq1)
        }
      })
      
      # Replace missing data with 0..
      af_vals[is.na(af_vals)] <- 0
      
      # Normalise for correlation matrix, otherwise keep as raw...
      if(normalise==TRUE){
        af_norm <- af_vals/sqrt(sum(af_vals^2))
      } else {
        af_norm <- af_vals
      }
      
      # If whole window is empty of snps, replace with 0...
      if(length(af_norm[is.na(af_norm)]) == length(af_norm)){
        af_norm <- rep(0,nrow(tmp_vcf@fix))
      }
      af_vec_0 <- data.frame(vec=af_norm)
      return(af_vec_0)
    })
    
    # Combine and export
    out <- matrix(nrow=length(vectors),ncol=nrow(vector_out[[1]]))
    for(i in 1:length(vectors)){
      out[i,] <- vector_out[[i]][,1]
    }
    colnames(out) <- sapply(rownames(vector_out[[1]]),function(snp){return(strsplit(snp,"_")[[1]][2])})
    rownames(out) <- names(vectors)
    return(out)
  },mc.cores=n_cores)
  
  # Fetch start and end pos as names
  names(window_list) <- sapply(1:length(winds),function(x){
    return(paste0(vcf@fix[,1][1],":",vcf@fix[,2][winds[x]],"-",vcf@fix[,2][winds2[x]]))
  })
  
  # Return
  return(window_list)
}

################################################################################
calc_linkage_matrices <- function(vcf=NULL, # Object of type vcfR, one chrom
                                  window_size=100, # window size in SNPs
                                  linkage_snps=NULL,
                                  n_cores=1){
  
  # Check for only one chr and make SNP windows
  if(length(unique(vcf@fix[,1])) > 1){
    stop("VCF contains more than one chromosome...")
  }
  winds <- seq(1,nrow(vcf@fix),window_size)
  winds2 <- winds+window_size-1
  winds2[length(winds2)] <- nrow(vcf@fix)
  
  # Trim any windows that are less than 10 SNPs...
  if(max(winds2)-max(winds) < 10){
    winds <- winds[1:(length(winds)-1)]
    winds2 <- winds2[1:(length(winds2)-1)]
  }
  
  # How many windows to estimate
  message(paste0("Calculating linkage matrices for ",length(winds)," windows with ",window_size," SNPs each"))
  
  # Run in parallel over all windows
  window_list <- mclapply(1:length(winds),function(x){
    
    # Subset the VCF and calculate linkage...
    tmp_vcf <- vcf[winds[x]:winds2[x],]
    genotypes <- t(extract.gt(tmp_vcf,return.alleles = T,convertNA = F))
    genotypes[genotypes == "./."] <- NA
    genotype_list <- lapply(1:ncol(genotypes),function(x){return(as.genotype(genotypes[,x]))})
    genotypes2 <- do.call(cbind.data.frame, genotype_list)
    
    # For now, we will limit to 50 SNPs per window...
    if(!(is.null(linkage_snps))){
      if(linkage_snps > window_snps){
        stop("Linkage SNP N is greater than window size")
      }
      ld_genotypes <- genetics::LD(genotypes2[,sample(1:ncol(genotypes2),linkage_snps)],)
    } else {
      ld_genotypes <- genetics::LD(genotypes2)
    }
    
    # Fetch window averages..
    ld_D <- median(na.omit(as.numeric(ld_genotypes$D)))
    ld_R <- median(na.omit(as.numeric(ld_genotypes$R^2)))
    return(c(ld_D,ld_R))
  },mc.cores = n_cores)
  
  
  # Fetch start and end pos as names
  names(window_list) <- sapply(1:length(winds),function(x){
    return(paste0(vcf@fix[,1][1],":",vcf@fix[,2][winds[x]],"-",vcf@fix[,2][winds2[x]]))
  })
  
  # Return
  return(window_list)
}

################################################################################
calc_null_AF_vectors <- function(vcf=NULL, # Object of type vcfR, one chrom
                                 window_size=100, # window size in SNPs
                                 popmap=NULL, # Popmap
                                 vectors=NULL, # List of all vectors from ancestor to descendent
                                 n_cores=1,
                                 perm=1000,
                                 normalise=TRUE){
  
  
  # Check for only one chr and make SNP windows
  if(length(unique(vcf@fix[,1])) > 1){
    stop("VCF contains more than one chromosome...")
  }
  
  # Draw windows randomly over whole chromosome
  max_wind <- nrow(vcf@fix)-window_size+1
  winds <- sort(sample(1:max_wind,perm))
  winds2 <- winds+window_size-1
  
  # How many windows to estimate
  message(paste0("Calculating Null AF vectors for ",perm," windows with ",window_size," SNPs each"))
  
  # Run in parallel over all windows
  window_list <- mclapply(1:length(winds),function(x){
    # print(x)
    
    # Shuffle the popmap for randomisations
    popmap[,1] <- sample(popmap[,1])
    
    # Make a "row" in the matrix for each vector
    vector_out <- lapply(vectors,function(vec){
      
      # Subset the VCF for each pop in the vector...
      tmp_vcf <- vcf[winds[x]:winds2[x],c(1,
                                          which(colnames(vcf@gt) %in% popmap[popmap[,2] == vec[1],1]),
                                          which(colnames(vcf@gt) %in% popmap[popmap[,2] == vec[2],1]))]
      
      # Calculate the allele frequencies and normalise
      af <- by(vcfR2loci(tmp_vcf),c(popmap[popmap[,2] == vec[1],2],
                                    popmap[popmap[,2] == vec[2],2]))
      
      # Get AF change
      af_vals <- sapply(af,function(i){
        # Mask bum windows
        if(ncol(i) == 0){
          return(0)
        } else {
          rownames(i) <- sort(vec)
          freq1 <- i[vec[1],1]/sum(i[vec[1],])
          freq2 <- i[vec[2],1]/sum(i[vec[2],])
          return(freq2 - freq1)
        }
      })
      
      # Replace missing data with 0..
      af_vals[is.na(af_vals)] <- 0
      
      # Normalise for correlation matrix, otherwise keep as raw...
      if(normalise==TRUE){
        af_norm <- af_vals/sqrt(sum(af_vals^2))
      } else {
        af_norm <- af_vals
      }      
      
      # If whole window is empty of snps, replace with 0...
      if(length(af_norm[is.na(af_norm)]) == length(af_norm)){
        af_norm <- rep(0,window_size)
      }
      af_vec_0 <- data.frame(vec=af_norm)
      return(af_vec_0)
    })
    
    # Combine and export
    out <- matrix(nrow=length(vectors),ncol=nrow(vector_out[[1]]))
    for(i in 1:length(vectors)){
      out[i,] <- vector_out[[i]][,1]
    }
    colnames(out) <- sapply(rownames(vector_out[[1]]),function(snp){return(strsplit(snp,"_")[[1]][2])})
    rownames(out) <- names(vectors)
    return(out)
  },mc.cores=n_cores)
  
  names(window_list) <- sapply(1:length(winds),function(x){
    return(paste0(vcf@fix[,1][1],":",vcf@fix[,2][winds[x]],"-",vcf@fix[,2][winds2[x]]))
  })
  
  # Return
  return(window_list)
}



################################################################################
# Calculate eigen information for each AF vector matrix...

eigen_analyse_vectors <- function(vector_input){
  C <- vector_input%*%t(vector_input)   ### Calculate matrix C ####
  eig <- eigen(C)  ### eigen decomposition of C ###
  vecs <- eig$vectors  ### eigenvectors (Q) of C ###
  vals <- eig$values
  a1 <- t(vecs)
  A <- t(vector_input) %*% solve(a1)
  
  # Tidy up
  rownames(vecs) <- rownames(vector_input)
  colnames(vecs) <- paste0("Eigenvector_",1:ncol(vecs))
  names(vals) <-  paste0("Eigenvector_",1:length(vals))
  
  # Return
  out <- list(vals,vecs,A)
  names(out) <- c("eigenvals","eigenvecs","A_matrix")
  return(out)
}

################################################################################
# Sum the eigenvalues...
# This function takes a list of all the eigen information and returns a list of summed values where each eigenvalue is the sum of itself the ones that preceed it...
sum_eigenvals <- function(eigen_res){
  
  # Count eigenvectors
  eigen_count <- length(eigen_res[[1]]$eigenvals)
  
  # Set up an out matrix
  out_mat <- matrix(nrow=length(eigen_res),ncol=eigen_count)
  rownames(out_mat) <- names(eigen_res)
  
  # Fill the mat
  for(i in 1:eigen_count){
    out_mat[,i] <- sapply(eigen_res,function(x){return(x$eigenvals[i])})
  }
  
  # Sum the mat
  for(i in 2:eigen_count){
    out_mat[,i] <- rowSums(out_mat[,c(i-1,i)])
  }
  
  # Clean
  colnames(out_mat) <- paste0("Eigenvalue_",1:eigen_count)
  
  # Return
  return(out_mat)
}

################################################################################
# Fetch the cutoff based on permutation...
find_null_cutoff <- function(null_res = NULL,
                             cutoffs = c(0.99)){
  # Do eigen analysis
  eigen_res <- lapply(null_res,eigen_analyse_vectors)
  
  # Count eigenvectors
  eigen_count <- length(eigen_res[[1]]$eigenvals)
  
  # Sum them
  null_sum <- sum_eigenvals(eigen_res)
  
  # Make an outmatrix...
  out_mat <- matrix(ncol=length(cutoffs),nrow=eigen_count)
  for(i in 1:ncol(out_mat)){
    for(j in 1:nrow(out_mat)){
      out_mat[j,i] <- quantile(null_sum[,j],probs = cutoffs[i])
    }
  }
  
  # Tidy and return
  colnames(out_mat) <- paste0(cutoffs*100,"%")
  rownames(out_mat) <- paste0("Eigenvector ",1:nrow(out_mat))
  return(out_mat)
}

################################################################################
# Look for significantly parallel windows
# parallel_eigenval_windows <- function(vector_input=NULL,
#                                       eigen_res=NULL,
#                                       window_size=100,
#                                       n_cores=1,
#                                       nsim=1000,
#                                       conf=0.999){
#   
#   # Estimate expected...
#   data_null <- matrix(nrow=nsim, ncol=nrow(vector_input[[1]]))
#   r<-rWishart::rWishart(nsim, window_size, diag(nrow(vector_input[[1]])), covariance = FALSE, simplify = "array")
#   for(i in 1:(nsim)){
#     C<-cov2cor(r[,,i])
#     eig<-eigen(C)
#     data_null[i,] <-eig$values
#   }
#   
#   # Now run over each window
#   res <- mclapply(1:length(vector_input),function(x){
#     # If last window, re-estimate the null due to different SNP N
#     if(x == length(vector_input)){
#       # Re-estimate expected...
#       data_null <- matrix(nrow=nsim, ncol=nrow(vector_input[[1]]))
#       r<-rWishart::rWishart(nsim, ncol(vector_input[[x]]), diag(nrow(vector_input[[1]])), covariance = FALSE, simplify = "array")
#       for(i in 1:(nsim)){
#         C<-cov2cor(r[,,i])
#         eig<-eigen(C)
#         data_null[i,] <-eig$values
#       }
#     } 
#     
#     # Otherwise analyse as normal
#     # How many signficant eigenvals?
#     signif <- sapply(1:length(eigen_res[[x]]$eigenvals),function(i){
#       ifelse(eigen_res[[x]]$eigenvals[i] > quantile(data_null[,i],probs=conf),return(TRUE),return(FALSE))
#     })
#     
#     # If none, return an empty cell
#     if(length(signif[signif == TRUE]) == 0){
#       return(NA)
#     } else {
#       return((1:length(eigen_res[[x]]$eigenvals))[signif])
#     }
#   })
#   return(res)
# }

# vector_input=AF_input
# eigen_res=eigen_res_test
# window_size=500
# n_cores=4
# nsim=1000
# conf=0.999


################################################################################
# Plot along a chr
eigenval_plot <- function(eigen_res=NULL,
                          cutoffs=NULL){
  
  # Count eigenvectors
  eigen_count <- length(eigen_res[[1]]$eigenvals)
  
  # Find window start and ends...
  pos_string <- sapply(names(eigen_res),function(str){return(strsplit(str,":")[[1]][2])})
  chr_string <- as.character(sapply(names(eigen_res),function(str){return(strsplit(str,":")[[1]][1])}))
  start <- as.integer(sapply(pos_string,function(str){return(strsplit(str,"-")[[1]][1])}))
  end <- as.integer(sapply(pos_string,function(str){return(strsplit(str,"-")[[1]][2])}))
  plot_dd <- data.frame(chr = chr_string,
                        pos = sort(c(start,end)))
  
  # Add in eigenvals
  eigenval_sums <- matrix(ncol=length(eigen_res[[1]]$eigenvals),nrow=nrow(plot_dd))
  for(i in 1:ncol(eigenval_sums)){
    eigenval_sums[,i] <- rep(as.numeric(sapply(eigen_res,function(j){return(j$eigenvals[i])})),
                             each=2)
  }
  colnames(eigenval_sums) <- paste0("Eigenvector ",1:ncol(eigenval_sums))
  
  # Merge and sum
  plot_dd <- cbind(plot_dd,eigenval_sums)
  for(i in 2:eigen_count){
    plot_dd[,paste0("Eigenvector ",i)] <- rowSums(plot_dd[,c(paste0("Eigenvector ",i-1),
                                                             paste0("Eigenvector ",i))])
  }
  
  # Visualise each eigenval along the chr
  plots <- lapply(colnames(eigenval_sums),function(val){
    
    p1 <- ggplot(plot_dd,aes(x=pos,y=plot_dd[,val]))+
      geom_line()+
      theme_bw()+
      labs(y=gsub("vector","value",val),x="Chr Position (Mb)")+
      ggtitle(val)+
      theme(axis.text=element_text(size=16),
            axis.title = element_text(size=18),
            title = element_text(size=20))+
      scale_x_continuous(breaks=seq(0,max(plot_dd$pos),2000000),
                         labels=seq(0,max(plot_dd$pos),2000000)/1000000)
    
    if(!(is.null(cutoffs))){
      p1 <- p1 + geom_hline(yintercept=cutoffs[val],colour="red2")
    }
    
    return(p1)
  })
  names(plots) <- colnames(plot_dd)[3:ncol(plot_dd)]
  return(plots)  
}

# eigen_res=eigen_res_test
# window_size=1000
# pos=pos_test
# nsim=1000
# conf=0.999

################################################################################
# Produce a PCA plot of suspected windows of interest
window_PCA_plot <- function(vcf = NULL,
                            window_id = NULL,
                            popmap = NULL){
  
  # Get the necessary information from the window_id
  chr <- strsplit(window_id,":")[[1]][1]
  bps <- strsplit(window_id,":")[[1]][2]
  start <- as.integer(strsplit(bps,"-")[[1]][1])
  end <- as.integer(strsplit(bps,"-")[[1]][2])
  
  # Subset the VCF
  sub_vcf <- vcf[vcf@fix[,1] == chr &
                   as.integer(vcf@fix[,2]) >= start &
                   as.integer(vcf@fix[,2]) <= end,]
  
  # Report
  message(paste0("Plotting PCA for window with ",nrow(sub_vcf@fix)," SNPs"))
  
  # Stop if we fail
  if(nrow(sub_vcf@fix)==0){
    stop("Error whilst subsetting, possible VCF does not contain chr or snps within the window")
  }
  
  # Change format
  pca_gen <- vcfR2genind(sub_vcf)
  x.gen<-tab(pca_gen,freq=TRUE,NA.method="mean")
  
  # Do PCA
  pca <- dudi.pca(x.gen,center=T,scale = F,scannf = FALSE, nf=10)
  scores<-pca$li
  scores$pop <- NA
  for(ind in popmap[,1]){
    scores[ind,"pop"] <- popmap[popmap[,1] == ind,2]
  }
  scores$pop_F <- factor(scores$pop,levels=unique(popmap$pop))
  
  return(s.class(pca$li,fac=scores$pop_F,col=funky(15)))
}

################################################################################
# Find significant windows
signif_eigen_windows <- function(eigen_res = NULL,
                                 cutoffs = NULL){
  
  if(is.numeric(cutoffs) == FALSE){
    stop("Error: cutoffs should be a numeric vector with a cutoff per eigenvector")
  }
  
  # Fetch sums
  eigen_sums <- sum_eigenvals(eigen_res)
  
  # Run for each eigenvector
  eigenvec_out <- lapply(1:ncol(eigen_sums),function(vec){
    
    # Fetch the significant
    signif_windows <- rownames(eigen_sums[eigen_sums[,vec] > cutoffs[vec],,drop=FALSE])
    
    # Return
    if(length(signif_windows) > 0){
      return(signif_windows)
    } else {
      return(NA)
    }
  })
  
  # Remove windows that are already significant based on previous eigenvalues...
  for(i in 2:length(eigenvec_out)){
    eigenvec_out[[i]] <- eigenvec_out[[i]][!(eigenvec_out[[i]] %in% unlist(lapply(1:(i-1),function(x){eigenvec_out[[x]]})))]
  }
  
  # Set names
  names(eigenvec_out) <- paste0("Eigenvector ",1:length(cutoffs))
  return(eigenvec_out)
}

################################################################################
# Summarise based on cutoffs
summarise_window_parallelism <- function(window_id,
                                         eigen_res,
                                         loading_cutoff = 0.3,
                                         eigenvector = 1){
  # Fetch loadings for eigenvectors
  loadings <- eigen_res[[window_id]]$eigenvecs[,eigenvector]
  
  # Output
  out <- data.frame(window_id=window_id,
                    parallel_lineages = length(loadings[abs(loadings) > loading_cutoff]))
  
  # Catch cases where all parallel but negative
  if(length(loadings[loadings < 0]) == length(loadings)){
    loadings <- -1 * loadings
  }
  
  # Now describe parallel/antiparallel
  out$parallel_pops <- paste(names(loadings[loadings >= loading_cutoff]),collapse = ",")
  out$antiparallel_pops <- paste(names(loadings[loadings <= (-1*loading_cutoff)]),collapse = ",")
  
  # Catch empty parallel_pops columns
  if(length(loadings[loadings >= loading_cutoff]) == 0){
    out$parallel_pops <- out$antiparallel_pops
    out$antiparallel_pops <- ""
  }
  
  return(out)
}

################################################################################
# Plot the eigenvalues for genome-wide/multiple chr
eigenval_plot_genome <- function(eigen_res_list=NULL,
                                 cutoffs=NULL){
  
  # Fetch chr to plot
  plot_chr <- unlist(lapply(eigen_res_list,function(x){return(names(x))}))
  plot_chr <- strsplit(plot_chr,":")
  plot_chr <- unique(sapply(plot_chr,function(x){return(x[1])}))
  
  # Fetch all the pos...
  total_plot_dd <- data.frame(rbindlist(lapply(plot_chr,function(chr){
    
    # Fetch eigen_res
    eigen_res <- eigen_res_list[[chr]]
    
    # Make plot_dd
    # Count eigenvectors
    eigen_count <- length(eigen_res[[1]]$eigenvals)
    
    # Find window start and ends...
    pos_string <- sapply(names(eigen_res),function(str){return(strsplit(str,":")[[1]][2])})
    chr_string <- as.character(sapply(names(eigen_res),function(str){return(strsplit(str,":")[[1]][1])}))
    start <- as.integer(sapply(pos_string,function(str){return(strsplit(str,"-")[[1]][1])}))
    end <- as.integer(sapply(pos_string,function(str){return(strsplit(str,"-")[[1]][2])}))
    plot_dd <- data.frame(chr = chr_string,
                          pos = sort(c(start,end)))
    
    # Add in eigenvals
    eigenval_sums <- matrix(ncol=length(eigen_res[[1]]$eigenvals),nrow=nrow(plot_dd))
    for(i in 1:ncol(eigenval_sums)){
      eigenval_sums[,i] <- rep(as.numeric(sapply(eigen_res,function(j){return(j$eigenvals[i])})),
                               each=2)
    }
    colnames(eigenval_sums) <- paste0("Eigenvector ",1:ncol(eigenval_sums))
    
    # Merge and sum
    plot_dd <- cbind(plot_dd,eigenval_sums)
    for(i in 2:eigen_count){
      plot_dd[,paste0("Eigenvector ",i)] <- rowSums(plot_dd[,c(paste0("Eigenvector ",i-1),
                                                               paste0("Eigenvector ",i))])
    }
    return(plot_dd)
  })))
  
  # Tidy up
  colnames(total_plot_dd)[3:ncol(total_plot_dd)] <- paste0("Eigenvector ",1:(ncol(total_plot_dd)-2))
  
  # Factorise according to input vector
  total_plot_dd$chr_F <- factor(total_plot_dd$chr,levels=chrs)
  
  
  # Visualise each eigenval along the chr
  plots <- lapply(grep("Eigen",colnames(total_plot_dd),value=T),function(val){
    
    p1 <- ggplot(total_plot_dd,aes(x=pos,y=total_plot_dd[,val]))+
      geom_line()+
      theme_minimal()+
      labs(y=gsub("vector","value",val),x="Genome Position")+
      ggtitle(val)+
      theme(axis.text=element_text(size=16),
            axis.title = element_text(size=18),
            title = element_text(size=20),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.spacing.x=unit(0.1, "lines"),
            strip.text = element_text(angle=90,hjust=1,size=14))+
      facet_wrap(~chr_F,nrow = 1,strip.position = "bottom",scales="free_x")
    
    if(!(is.null(cutoffs))){
      p1 <- p1 + geom_hline(yintercept=cutoffs[val],colour="red2")
    }
  })
  names(plots) <- grep("Eigen",colnames(total_plot_dd),value=T)
  return(plots)  
}



################################################################################
# haplotype_parallel <- function(vcf=NULL,
#                                windows=NULL,
#                                popmap=NULL,
#                                ancestral=NULL,
#                                n_cores=1){
#   
#   # Define our focal pops
#   focals <- grep(ancestral,unique(popmap[,2]),invert=T,value=T)
#   
#   library(vcfR)
#   library(adegenet)
#   library(dplyr)
#   library(Morpho)
#   library(reshape2)
#   library(parallel)
#   
#   window_out <- mclapply(1:length(windows),function(N){
#     
#     # Print start
#     print(paste0("Starting window ",windows[N]))
#     
#     # Read in the window
#     system(paste0("bcftools view -r ",windows[N]," ",vcf," > ",N,".vcf"))
#     wind_vcf <- read.vcfR(paste0(N,".vcf"),verbose = F)
#     wind_vcf <- wind_vcf <- wind_vcf[,c(1,which(colnames(wind_vcf@gt) %in% popmap[,1]))]
#     system(paste0("rm -f ",N,".vcf"))
#     
#     # Convert
#     wind_genind <- vcfR2genind(wind_vcf)
#     x.gen<-tab(wind_genind,freq=TRUE,NA.method="mean")
#     
#     # Do PCA
#     pca <- dudi.pca(x.gen,center=T,scale = F,scannf = FALSE, nf=10)
#     scores <- pca$li
#     #s.label(pca$li)
#     
#     # Average according to popmap
#     scores <- scores[popmap[,1],]
#     scores$pop <- popmap[,2]
#     avg_scores <- data.frame(scores %>% group_by(pop) %>%
#                                summarise_all(mean))
#     
#     # Scale vectors according to weighting...
#     for(i in 2:ncol(avg_scores)){
#       avg_scores[,i] <- avg_scores[,i] * (pca$eig/sum(pca$eig))[i]
#     }
#     
#     # Fetch all the vectors
#     vectors <- lapply(focals,function(x){
#       vec1 <- as.vector(t(avg_scores[avg_scores$pop == ancestral,2:ncol(avg_scores)] - avg_scores[avg_scores$pop == x,2:ncol(avg_scores)]))
#     })
#     
#     # Combinations
#     to_compare <- combn(1:length(focals),2)
#     
#     # Compare the vectors
#     theta_mat <- matrix(nrow=length(focals),ncol=length(focals))
#     length_mat <- matrix(nrow=length(focals),ncol=length(focals))
#     
#     # Fill the matrix
#     for(i in 1:ncol(to_compare)){
#       comp1 <- to_compare[1,i]
#       comp2 <- to_compare[2,i]
#       
#       # Calc angle
#       angle.deg <- angle.calc(vectors[[comp1]],vectors[[comp2]])*(180/pi)
#       
#       # Calc difference in length
#       get.vectlength <- function(vect){
#         return(sqrt(sum(vect^2, na.rm = TRUE)))
#       }
#       deltaL<-abs(round(get.vectlength(vectors[[comp1]]) - get.vectlength(vectors[[comp2]]), 3))
#       
#       # Fill the matrices
#       theta_mat[comp1,comp2] <- angle.deg
#       length_mat[comp1,comp2] <- deltaL
#     }
#     
#     # Tidy up
#     colnames(theta_mat) <- focals
#     rownames(theta_mat) <- focals
#     colnames(length_mat) <- focals
#     rownames(length_mat) <- focals
#     
#     # Make alternative outputs
#     theta_out <- na.omit(melt(theta_mat))
#     colnames(theta_out) <- c("pop1","pop2","theta")
#     length_out <- na.omit(melt(length_mat))
#     colnames(length_out) <- c("pop1","pop2","deltaL")
#     theta_out$window <- windows[N]
#     length_out$window <- windows[N]
#     
#     # Return everything on a per window basis
#     return(list(theta_mat,length_mat,
#                 theta_out,length_out))
#     
#   },mc.cores=n_cores)
#   return(window_out)
# }
# 
# # Function for fetching
# fetch_parallel_output <- function(parallel_res=NULL,
#                                   stat="theta"){
#   library(data.table)
#   if(stat == "theta"){
#     data.frame(rbindlist(lapply(parallel_res,function(x){return(x[[3]])})))
#   } else if (stat == "deltaL"){
#     data.frame(rbindlist(lapply(parallel_res,function(x){return(x[[4]])})))
#   } else {
#     stop("Not a valid option!!")
#   }
# }
# 
# # # Run the function
# # parallel_stats <- haplotype_parallel(vcf = vcf_test,
# #                                      windows = test_winds,
# #                                      popmap = test_popmap,
# #                                      ancestral ="GH")
# # 
# # # Get results
# # theta_res <- fetch_parallel_output(parallel_stats,stat = "theta")
# # deltaL_res <- fetch_parallel_output(parallel_stats,stat = "deltaL")
# # 
# # # Plot the comparisons
# # library(ggplot2)
# # theta_res$comp <- paste0(theta_res$pop1,"_",theta_res$pop2)

################################################################################
#### pop_structure_AF ####
# This function attempts to recreate population structure by permuting a set of random snps equal in number to window_snps and creating an eigenvalue distribution to test against...
pop_structure_AF <- function(vcf,
                             vectors,
                             popmap,
                             window_snps,
                             permutations=100,
                             normalise=TRUE,
                             n_cores=1){
  
  # Permute over N permutations
  permute_out <- mclapply(1:permutations,function(perm){
    
    # Downsample vcf randomly
    perm_vcf <- vcf[sort(sample(1:nrow(vcf),window_snps)),]
    
    # Make a "row" in the matrix for each vector
    vector_out <- lapply(vectors,function(vec){
      
      # Subset the VCF for each pop in the vector...
      tmp_vcf <- perm_vcf[,c(1,which(colnames(vcf@gt) %in% popmap[popmap[,2] == vec[1],1]),
                             which(colnames(vcf@gt) %in% popmap[popmap[,2] == vec[2],1]))]
      
      # Calculate the allele frequencies and normalise
      af <- by(vcfR2loci(tmp_vcf),c(popmap[popmap[,2] == vec[1],2],
                                    popmap[popmap[,2] == vec[2],2]))
      
      # Get AF change
      af_vals <- sapply(af,function(i){
        freq1 <- i[vec[1],1]/sum(i[vec[1],])
        freq2 <- i[vec[2],1]/sum(i[vec[2],])
        return(freq2 - freq1)
      })
      
      # Replace missing data with 0..
      af_vals[is.na(af_vals)] <- 0
      
      # Normalise for correlation matrix, otherwise keep as raw...
      if(normalise==TRUE){
        af_norm <- af_vals/sqrt(sum(af_vals^2))
      } else {
        af_norm <- af_vals
      }
      
      # If whole window is empty of snps, replace with 0...
      if(length(af_norm[is.na(af_norm)]) == length(af_norm)){
        af_norm <- rep(0,window_size)
      }
      af_vec_0 <- data.frame(vec=af_norm)
      return(af_vec_0)
    })
    
    # Combine and export
    out <- matrix(nrow=length(vectors),ncol=nrow(vector_out[[1]]))
    for(i in 1:length(vectors)){
      out[i,] <- vector_out[[i]][,1]
    }
    colnames(out) <- sapply(rownames(vector_out[[1]]),function(snp){return(strsplit(snp,"_")[[1]][2])})
    rownames(out) <- names(vectors)
    return(out)
  },mc.cores=n_cores)
  
  # Name permuted AF_matrices
  names(permute_out) <- paste0("perm_",1:permutations)
  return(permute_out) 
}

# # Test
# vcf <- vcf_in
# window_snps = 500
# permutations = 500
# n_cores=6
# vectors = test_vectors
# popmap = test_popmap
# 
# eigen_null <- lapply(permute_out,function(af){
#   eigen_analyse_vectors(af,normalise = TRUE)
# })
# 
# hist(sapply(eigen_null,function(x){return(x$eigenvals[1])}))

################################################################################
#### mds_transform_AF ####
# This script uses functions from lostruct packages to MDS transform our AF eigen matrices...
mds_transform_AF <- function(vector_input,
                             n_cores,
                             npc=length(vector_input[[1]]$eigenvals)){
  
  # First check for lostruct
  if(!("lostruct" %in% rownames(installed.packages()))){
    stop("R package lostruct is required but not installed, visit https://github.com/petrelharp/local_pca")
  } 
  library(lostruct)
  
  # Take the eigen windows and transform prior to mds scaling
  # Format needs to be single matrix of row per window, col1 = sum of squares, col2 onwards is eigenvector loadings followed by values
  k = nrow(vector_input[[1]])
  eigen.reformat <- matrix(ncol=(1+k+k^2),nrow=length(vector_input))
  for(i in 1:nrow(eigen.reformat)){
    x <- vector_input[[i]]
    tmp1 <- sum(as.vector(x^2))
    tmp2 <- eigen_analyse_vectors(x)$eigenvals
    tmp3 <- eigen_analyse_vectors(x)$eigenvecs
    eigen.reformat[i,] <- c(tmp1,tmp2,tmp3)
  }
  
  # Compute distance matrix
  attr(eigen.reformat,"npc") <- k
  mds_out <- lostruct::pc_dist(eigen.reformat,mc.cores = n_cores)
  
  # Return
  return(mds_out)
}

# # Test
# 
# windist <- mds_transform_AF(AF_input,6)
# fit2d <- cmdscale(windist, eig=TRUE, k=3)
# to_plot <- data.frame(fit2d$points)
# colnames(to_plot) <- paste0("MDS",1:ncol(to_plot))
# BPs <- names(AF_input)
# BPs <- strsplit(BPs,":")
# BPs <- unlist(lapply(BPs,function(x){return(x[2])}))
# BPs <- strsplit(BPs,"-")
# BPs <- unlist(lapply(BPs,function(x){return(mean(as.integer(x)))}))
# 
# to_plot$BP <- BPs[1:nrow(to_plot)]
# 
# # Plot scatter
# ggplot(to_plot,aes(MDS1,MDS2,colour=BP))+geom_point()
# 
# # Plot MDS separately
# mds_plot <- lapply(1:3,function(mds){
# ggplot(to_plot,aes(BP,y=to_plot[,mds]))+
#     geom_point()+
#     ylab(paste0("MDS",mds))+
#     xlab("BP")
# })
# cowplot::plot_grid(plotlist = mds_plot,
#           ncol=1)