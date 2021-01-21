# this script calculates various allele frequency 'stats'
library(data.table)
### mean/median/min-max of absolute AF changes
gh<-fread('output/allele_freq/GHP_minor.txt',header=T)
ll<-fread('output/allele_freq/LL_minor.txt',header=T)
ul<-fread('output/allele_freq/UL_minor.txt',header=T)
c<-fread('output/allele_freq/C_minor.txt',header=T)
t<-fread('output/allele_freq/T_minor.txt',header=T)
gl<-fread('output/allele_freq/GLP_minor.txt',header=T)

gh$ID<-paste0(gh$CHROM,':',gh$POS)
ll$ID<-paste0(ll$CHROM,':',ll$POS)
ul$ID<-paste0(ul$CHROM,':',ul$POS)
c$ID<-paste0(c$CHROM,':',c$POS)
t$ID<-paste0(t$CHROM,':',t$POS)
gl$ID<-paste0(gl$CHROM,':',gl$POS)

# create a df for each HP/LP pair
ghll2<-data.frame(LL_MAF=ll$minor_freq,GH_MAF=gh$GHP_minor_FRQ, ID=ll$ID)
ghul2<-data.frame(UL_MAF=ul$minor_freq,GH_MAF=gh$GHP_minor_FRQ, ID=ul$ID)
ghc2<-data.frame(C_MAF=c$minor_freq,GH_MAF=gh$GHP_minor_FRQ, ID=c$ID)
ght2<-data.frame(T_MAF=t$minor_freq,GH_MAF=gh$GHP_minor_FRQ, ID=t$ID)
ghgl2<-data.frame(GL_MAF=gl$minor_freq,GH_MAF=gh$GHP_minor_FRQ, ID=gl$ID)

# calculate absolute allele frequency changes between each LP and GHP
ghll2$absDAF<-abs(ghll2$LL_MAF-ghll2$GH_MAF)
ghul2$absDAF<-abs(ghul2$UL_MAF-ghul2$GH_MAF)
ghc2$absDAF<-abs(ghc2$C_MAF-ghc2$GH_MAF)
ght2$absDAF<-abs(ght2$T_MAF-ght2$GH_MAF)
ghgl2$absDAF<-abs(ghgl2$GL_MAF-ghgl2$GH_MAF)

# create empty df
af_mat<-data.frame(matrix(ncol=7,nrow = 5))  # nrow = the nr of pairwise comparisons
colnames(af_mat)<-c("Population","Fixed GHP/LP","GHP MAF/fixed LP","Unchanged MAF","mean |dAF|","median |dAF|","min-max")

# start population the df
af_mat$Population<-c("ILL","IUL","IC","IT","GLP")

# fill in the mean and median absolute AF changes
af_mat$`mean |dAF|`<-c(round(mean(ghll2$absDAF),digits=3),round(mean(ghul2$absDAF),digits=3),round(mean(ghc2$absDAF),digits=3),
                       round(mean(ght2$absDAF),digits=3),round(mean(ghgl2$absDAF),digits=3))
af_mat$`median |dAF|`<-c(round(median(ghll2$absDAF),digits=3),round(median(ghul2$absDAF),digits=3),round(median(ghc2$absDAF),digits=3),
                         round(median(ght2$absDAF),digits=3),round(median(ghgl2$absDAF),digits=3))

# fill in the min & max
ll_min<-paste0(min(ghll2$absDAF),'-',round(max(ghll2$absDAF),digits=3))
ul_min<-paste0(min(ghul2$absDAF),'-',round(max(ghul2$absDAF),digits=3))
c_min<-paste0(min(ghc2$absDAF),'-',round(max(ghc2$absDAF),digits=3))
t_min<-paste0(min(ght2$absDAF),'-',round(max(ght2$absDAF),digits=3))
gl_min<-paste0(min(ghgl2$absDAF),'-',round(max(ghgl2$absDAF),digits=3))

af_mat$`min-max`<-c(ll_min,ul_min,c_min,t_min,gl_min)

# calculate unchanged AFs
ll_same<-ghll2[ghll2$LL_MAF==ghll2$GH_MAF ,]
ul_same<-ghul2[ghul2$UL_MAF==ghul2$GH_MAF ,]
c_same<-ghc2[ghc2$C_MAF==ghc2$GH_MAF ,]
t_same<-ght2[ght2$T_MAF==ght2$GH_MAF ,]
gl_same<-ghgl2[ghgl2$GL_MAF==ghgl2$GH_MAF ,]

af_mat$`Unchanged MAF`<-c(nrow(ll_same),nrow(ul_same),nrow(c_same),nrow(t_same),nrow(gl_same))

### Find the nr of fixed differences between GHP and each LP
# read in the GHP file
gh_maf<-fread('data/allele_freq/GHP_AF.txt',header=T)
gh_maf$ID<-paste0(gh_maf$CHROM,':',gh_maf$POS)
gh_maf<-gh_maf[,c(1,2,5,6,7,8,9)]
colnames(gh_maf)<-c('CHROM','POS','GH_REF','GH_REF_FRQ','GH_ALT','GH_ALT_FRQ','ID')

# which sites are fixed in GHP
gh_fix_ref<-gh_maf[gh_maf$GH_REF_FRQ==1,]
gh_fix_alt<-gh_maf[gh_maf$GH_ALT_FRQ==1,]
#gh_fix2<-gh_maf[gh_maf$GH_ALT_FRQ==1,]

# read in the GLP data
gl_maf<-fread('data/allele_freq/GLP_AF.txt',header=T)
gl_maf$ID<-paste0(gl_maf$CHROM,':',gl_maf$POS)
gl_maf<-gl_maf[order(gl_maf$CHROM,decreasing = F),]

# Find how many sites were the minor allele in GH and are now fixed in LP
gl_maf<-gl_maf[order(gl_maf$CHROM,decreasing = F),]
gl1<-gl_maf[gl_maf$REF==gh$minor_allele,]
gl1b<-gl1[gl1$REF_FRQ==1,]
gl2<-gl_maf[gl_maf$ALT==gh$minor_allele,]
gl2b<-gl2[gl2$ALT_FRQ==1,]
gl_min_fix<-nrow(gl1b)+nrow(gl2b)

# fixed between HP & LP
gl_ref<-gl_maf[gl_maf$ID%in%gh_fix_ref$ID,]
gl_ref<-gl_ref[gl_ref$REF_FRQ==0,]
gl_ref2<-gl_maf[gl_maf$ID%in%gh_fix_alt$ID,]
gl_ref2<-gl_ref2[gl_ref2$ALT_FRQ==0,]
gl_fixed<-nrow(gl_ref)+nrow(gl_ref2)

# Read in LL data
ll_maf<-fread('data/allele_freq/ILL_AF.txt',header=T)
ll_maf$ID<-paste0(ll_maf$CHROM,':',ll_maf$POS)

# minor in GH/fixed in LP
ll_maf<-ll_maf[order(ll_maf$CHROM,decreasing = F),]
ll1<-ll_maf[ll_maf$REF==gh$minor_allele,]
ll1b<-ll1[ll1$REF_FRQ==1,]
ll2<-ll_maf[ll_maf$ALT==gh$minor_allele,]
ll2b<-ll2[ll2$ALT_FRQ==1,]
#ll1c<-ll1[ll1$REF_FRQ==0,]
ll_min_fix<-nrow(ll1b)+nrow(ll2b)

# fixed between HP & LP
ll_ref<-ll_maf[ll_maf$ID%in%gh_fix_ref$ID,]
ll_ref<-ll_ref[ll_ref$REF_FRQ==0,]
ll_ref2<-ll_maf[ll_maf$ID%in%gh_fix_alt$ID,]
ll_ref2<-ll_ref2[ll_ref2$ALT_FRQ==0,]
ll_fixed<-nrow(ll_ref)+nrow(ll_ref2)

# UL
ul_maf<-fread('data/allele_freq/IUL_AF.txt',header=T)
ul_maf$ID<-paste0(ul_maf$CHROM,':',ul_maf$POS)

# minor in GH/fixed in LP
ul_maf<-ul_maf[order(ul_maf$CHROM,decreasing = F),]
ul1<-ul_maf[ul_maf$REF==gh$minor_allele,]
ul1b<-ul1[ul1$REF_FRQ==1,]
ul2<-ul_maf[ul_maf$ALT==gh$minor_allele,]
ul2b<-ul2[ul2$ALT_FRQ==1,]
#ul1c<-ul1[ul1$REF_FRQ==0,]
ul_min_fix<-nrow(ul1b)+nrow(ul2b)

# fixed between HP & LP
ul_ref<-ul_maf[ul_maf$ID%in%gh_fix_ref$ID,]
ul_ref<-ul_ref[ul_ref$REF_FRQ==0,]
ul_ref2<-ul_maf[ul_maf$ID%in%gh_fix_alt$ID,]
ul_ref2<-ul_ref2[ul_ref2$ALT_FRQ==0,]
ul_fixed<-nrow(ul_ref)+nrow(ul_ref2)

# C
c_maf<-fread('data/allele_freq/IC_AF.txt',header=T)
c_maf$ID<-paste0(c_maf$CHROM,':',c_maf$POS)

# minor in GH/fixed in LP
c_maf<-c_maf[order(c_maf$CHROM,decreasing = F),]
c1<-c_maf[c_maf$REF==gh$minor_allele,]
c1b<-c1[c1$REF_FRQ==1,]
c2<-c_maf[c_maf$ALT==gh$minor_allele,]
c2b<-c2[c2$ALT_FRQ==1,]
#c1c<-c1[c1$REF_FRQ==0,]
c_min_fix<-nrow(c1b)+nrow(c2b)

# fixed between HP & LP
c_ref<-c_maf[c_maf$ID%in%gh_fix_ref$ID,]
c_ref<-c_ref[c_ref$REF_FRQ==0,]
c_ref2<-c_maf[c_maf$ID%in%gh_fix_alt$ID,]
c_ref2<-c_ref2[c_ref2$ALT_FRQ==0,]
c_fixed<-nrow(c_ref)+nrow(c_ref2)

# T
t_maf<-fread('data/allele_freq/IT_AF.txt',header=T)
t_maf$ID<-paste0(t_maf$CHROM,':',t_maf$POS)

# minor in GH/fixed in LP
t_maf<-t_maf[order(t_maf$CHROM,decreasing = F),]
t1<-t_maf[t_maf$REF==gh$minor_allele,]
t1b<-t1[t1$REF_FRQ==1,]
t2<-t_maf[t_maf$ALT==gh$minor_allele,]
t2b<-t2[t2$ALT_FRQ==1,]
#t1c<-t1[t1$REF_FRQ==0,]
t_min_fix<-nrow(t1b)+nrow(t2b)

# fixed between HP & LP
t_ref<-t_maf[t_maf$ID%in%gh_fix_ref$ID,]
t_ref<-t_ref[t_ref$REF_FRQ==0,]
t_ref2<-t_maf[t_maf$ID%in%gh_fix_alt$ID,]
t_ref2<-t_ref2[t_ref2$ALT_FRQ==0,]
t_fixed<-nrow(t_ref)+nrow(t_ref2)

# populate the rest of the table
af_mat$`GHP MAF/fixed LP`<-c(ll_min_fix,ul_min_fix,c_min_fix,t_min_fix,gl_min_fix)
af_mat$`Fixed GHP/LP`<-c(ll_fixed,ul_fixed,c_fixed,t_fixed,gl_fixed)


write.table(af_mat, "output/allele_freq/AF_fixed_table.txt", quote = F, sep = '\t', row.names = F)
