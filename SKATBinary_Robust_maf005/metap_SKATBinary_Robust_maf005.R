# Meta-analysis (p-values obtained via SKATBinary_Robust)

rm(list = ls())
library(dplyr)
library(tidyr)
#install.packages("metap")
library(metap)


#######################################################################################################################
#  ALL

filepath.out <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/SKATO_test_results_SKATBinary_Robust_BG_maf005/"
pop1 <- "ALLcasecontrol_cohort1"; pop2 <- "ALLcasecontrol_cohort2"
c1 <- read.table(paste0(filepath.out, pop1, "/summary_results/allResults.txt"),header = TRUE)
c2 <- read.table(paste0(filepath.out, pop2, "/summary_results/allResults.txt"),header = TRUE)
tab.set.c1 <- table(c1$SetID)
dup.setID.c1 <- rownames(tab.set.c1[tab.set.c1 > 1])
tmp.pval <- c()
for (s in dup.setID.c1) {
  print( c1$Pvalue[c1$SetID == s])
  tmp.pval <- c(tmp.pval, unlist(c1$Pvalue[c1$SetID == s]))
}
min(tmp.pval)
c1 <- c1[!(c1$SetID %in% dup.setID.c1), ]
combined.Pvalues <- c1 %>% left_join(c2, by=c("SetID"))

n.meta <- dim(combined.Pvalues)[1]
print(paste0("There are ", n.meta, " effective tests in Meta, and the corresponding cutoff value is ", 0.05/n.meta))

# sumz
# By default the weights are equal. In the absence of effect sizes (in which case a method for combining effect
# sizes would be more appropriate anyway) best results are believed to be obtained with weights
# proportional to the square root of the sample sizes (Zaykin 2011)
n1 <- 2425; n2 <- 555
weights <- c(sqrt(n1), sqrt(n2))
p.value.sumz <- rep(NA, n.meta)
for (i in 1:n.meta) {
  pvalues <- unlist(combined.Pvalues[i, c(2,3)])
  if (is.na(pvalues[2])) {
    p.value.sumz[i] <- pvalues[1]
    next;
  }
  
  temp <- try(sumz(pvalues, weights = weights))
  if (is.na(temp$p)) {
    p.value.sumz[i] <- pvalues[1]
  } else {
    p.value.sumz[i] <- temp$p
  }
}
#--------------------------------------------------------------------------------------------#
res <- data.frame("SetID"=combined.Pvalues$SetID, "meta_sumz_Pvalue"=p.value.sumz)

res[res$SetID == "EVC2",]
res[res$meta_sumz_Pvalue < 0.05 / 12741, ]



###############################################################################################
# qqplot 
library(qqman)

# ALL
## cohort1
cutoff.ALL <- 0.05 / 12741; transformed.cutoff.ALL <- - log10(cutoff.ALL)

## sumz
setwd("/user/xinyugao/SKATO_Project/SKATBinary_Robust_maf005/")
png(paste0("ALL.meta_sumz.SKATBinary_Robust_maf005.png"))
qq(res$meta_sumz_Pvalue, main=paste("QQ-plot of ALL Meta Pvalues via Stouffer's with cutoff val=", signif(cutoff.ALL, 2)), 
   cex.main=1, ylim =c(0,-log10(1e-07) ))
abline(h=-log10(cutoff.ALL), col="blue")
dev.off()


## Manhattan plots
filepath <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/SKATO_test_results_SKATBinary/Chr_BP_Genes/"
chr_bp.info <- read.table(paste0(filepath, "chr_BP_Genes.txt"), header = TRUE)
colnames(chr_bp.info)[1] <-"SetID"
head(chr_bp.info)


res.dat <- res %>% left_join(chr_bp.info, by=c("SetID"))
res.dat$Chr <- as.numeric(res.dat$Chr)
setwd("/user/xinyugao/SKATO_Project/SKATBinary_Robust_maf005/")
png(paste0("ALL.meta_sumz.manhattan_SKATBinary_Robust_maf005.png"), units="px", height=8000, width=15000, res=1015)
manhattan(res.dat, chr='Chr', bp='BP', p="meta_sumz_Pvalue", snp='SetID', main="Manhattan plot of Meta Pvalues for ALL",
          col=c("dodgerblue2","darkblue"), cex.axis=.7, las=2, 
          annotatePval = cutoff.ALL, annotateTop = F, 
          genomewideline = F, suggestiveline=F)
abline(h=-log10(cutoff.ALL), col="red")
dev.off()





# cohort1
res.dat.cohort1 <- c1 %>% left_join(chr_bp.info, by = c("SetID"))
res.dat.cohort1[which.minn(res.dat.cohort1$Pvalue, 10),]
# cohort2
res.dat.cohort2 <- c2 %>% left_join(chr_bp.info, by = c("SetID"))
res.dat.cohort2[which.minn(res.dat.cohort2$Pvalue, 10),]
# meta
res.dat[which.minn(res.dat$meta_sumz_Pvalue, 10),]






filepath.cohort1 <- paste0("/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/SKAT_META/META_SSD_RareVariants_maf005/ALLcasecontrol_cohort1")
filepath.cohort2 <- paste0("/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/SKAT_META/META_SSD_RareVariants_maf005/ALLcasecontrol_cohort2")
File.MSSD.vec <- c(paste0(filepath.cohort1, "/File.MSSD"),
                   paste0(filepath.cohort2, "/File.MSSD"))
File.MInfo.vec <- c(paste0(filepath.cohort1, "/File.MInfo"),
                    paste0(filepath.cohort2, "/File.MInfo"))
res.open_MSSD <- Open_MSSD_File_2Read(File.MSSD.vec = File.MSSD.vec, 
                                      File.MInfo.vec = File.MInfo.vec)

Info_MAF_Missing_Score_Allele12.ALL.cohort1  <- res.open_MSSD$EachInfo[[1]]$Info
Info_MAF_Missing_Score_Allele12.ALL.cohort2  <- res.open_MSSD$EachInfo[[2]]$Info
EVC2.cohort1 <- Info_MAF_Missing_Score_Allele12.ALL.cohort1[Info_MAF_Missing_Score_Allele12.ALL.cohort1$SetID=="EVC2",]
EVC2.cohort2 <- Info_MAF_Missing_Score_Allele12.ALL.cohort2[Info_MAF_Missing_Score_Allele12.ALL.cohort2$SetID=="EVC2",]
dim(EVC2.cohort1)
EVC2.cohort1
EVC2.cohort2


signif(EVC2.cohort1$Score, 3)

cbind(signif(EVC2.cohort2$Score, 3), signif(EVC2.cohort2$MAF, 2))


