rm(list = ls())
library(dplyr)
library(tidyr)
#install.packages("metap")
library(metap)


#######################################################################################################################
#  ALL

filepath.out <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/SKATO_test_results_SKATBinary_Robust_BG_maf005/"
pop1 <- "ALLcasecontrol_cohort1"; pop2 <- "ALLcasecontrol_cohort2"
Tcell.c1 <- read.table(paste0(filepath.out, "Tcell/",pop1, "/summary_results/allResults.txt"),header = TRUE)
dim(Tcell.c1) # 12490
Tcell.c2 <- read.table(paste0(filepath.out, "Tcell/",pop2, "/summary_results/allResults.txt"),header = TRUE)
dim(Tcell.c2) # 9525 
Bcell.c1 <- read.table(paste0(filepath.out, "Bcell/",pop1, "/summary_results/allResults.txt"),header = TRUE)
dim(Bcell.c1) # 12639
Bcell.c2 <- read.table(paste0(filepath.out, "Bcell/",pop2, "/summary_results/allResults.txt"),header = TRUE)
dim(Bcell.c2)  # 9678


check_dup_setID <- function(c1) {
  tab.set.c1 <- table(c1$SetID)
  dup.setID.c1 <- rownames(tab.set.c1[tab.set.c1 > 1])
  tmp.pval <- c()
  for (s in dup.setID.c1) {
    #print( c1$Pvalue[c1$SetID == s])
    tmp.pval <- c(tmp.pval, unlist(c1$Pvalue[c1$SetID == s]))
    
  }
  
  print(min(tmp.pval, na.rm = TRUE))
  c1 <- c1[!(c1$SetID %in% dup.setID.c1), ]
  return(c1)
}
Tcell.c1 <- check_dup_setID(Tcell.c1)
Bcell.c1 <- check_dup_setID(Bcell.c1)
dim(Tcell.c1) # 12389



get_meta_pvals <- function(n1, n2, combined.Pvalues) {
  weights <- c(sqrt(n1), sqrt(n2))
  n.meta <- dim(combined.Pvalues)[1]
  p.value.sumz <- rep(NA, n.meta)
  for (i in 1:n.meta) {
    pvalues <- unlist(combined.Pvalues[i, 2:3])
    if (is.na(pvalues[2])) {
      p.value.sumz[i] <- as.numeric(pvalues[1])
      next;
    }
    temp <- try(sumz(pvalues, weights = weights))
    if (is.na(temp$p)) {
      p.value.sumz[i] <- as.numeric(pvalues[1])
    } else {p.value.sumz[i] <- temp$p}
  }
  return(p.value.sumz)
}



n1.Tcell <- 2046; n2.Tcell <- 518
combined.Pvalues.Tcell <- Tcell.c1 %>% left_join(Tcell.c2, by=c("SetID"))
metap.Tcell <- get_meta_pvals(n1 = n1.Tcell, n2 = n2.Tcell, combined.Pvalues = combined.Pvalues.Tcell)
res.meta.Tcell <- data.frame("Gene"=combined.Pvalues.Tcell$SetID, "meta_sumz_Pvalue"=metap.Tcell)


n1.Bcell <- 2293; n2.Bcell <- 549
combined.Pvalues.Bcell <- Bcell.c1 %>% left_join(Bcell.c2, by=c("SetID"))
metap.Bcell <- get_meta_pvals(n1 = n1.Bcell, n2 = n2.Bcell, combined.Pvalues = combined.Pvalues.Bcell)
res.meta.Bcell <- data.frame("Gene"=combined.Pvalues.Bcell$SetID, "meta_sumz_Pvalue"=metap.Bcell)



library(qqman)

n.Tcell.eff1 <- 12490; n.Tcell.eff2 <- 9525 
n.Bcell.eff1 <- 12639; n.Bcell.eff2 <- 9678
cutoff.Tcell <- 0.05/n.Tcell.eff1; cutoff2.Tcell <- 0.05/n.Tcell.eff2
cutoff.Bcell <- 0.05/n.Bcell.eff1; cutoff2.Bcell <- 0.05/n.Bcell.eff2

# QQ-plot
# T-cell
setwd("/user/xinyugao/SKATO_Project/SKATBinary_Robust_maf005/Tcell/")
png(paste0("Tcell.meta_sumz.SKATBinary_Robust_maf005.png")) 
qq(res.meta.Tcell$meta_sumz_Pvalue, main=paste("QQ-plot of T-cell Meta Pvalues with cutoff=", signif(cutoff.Tcell, 2)), cex.main=1, 
   ylim =c(0,-log10(1e-06))  )
abline(h=-log10(cutoff.Tcell), col="blue")
dev.off()

# B-cell
setwd("/user/xinyugao/SKATO_Project/SKATBinary_Robust_maf005/Bcell/")
png(paste0("Bcell.meta_sumz.SKATBinary_Robust_maf005.png")) 
qq(res.meta.Bcell$meta_sumz_Pvalue, main=paste("QQ-plot of B-cell Meta Pvalues with cutoff=", signif(cutoff.Bcell, 2)), cex.main=1, 
   ylim =c(0,-log10(1e-06))  )
abline(h=-log10(cutoff.Bcell), col="blue")
dev.off()


# Manhattan plot
filepath <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/SKATO_test_results_SKATBinary/Chr_BP_Genes/"
chr_bp.info <- read.table(paste0(filepath, "chr_BP_Genes.txt"), header = TRUE)
colnames(chr_bp.info)[1] <-"SetID"
head(chr_bp.info)

# T-cell
res.dat.Tcell <- res.meta.Tcell %>% left_join(chr_bp.info, by=c("Gene" = "SetID"))
res.dat.Tcell$Chr <- as.numeric(res.dat.Tcell$Chr)
setwd("/user/xinyugao/SKATO_Project/SKATBinary_Robust_maf005/Tcell/")
png(paste0("Tcell.meta_sumz.manhattan_SKATBinary_Robust_maf005.png"), units="px", height=8000, width=15000, res=1015)
manhattan(res.dat.Tcell, chr='Chr', bp='BP', p="meta_sumz_Pvalue", snp='Gene', main="Manhattan plot of Meta Pvalues for T-cell",
          col=c("dodgerblue2","darkblue"), cex.axis=.7, las=2, ylim=c(0, 6),
          annotatePval = cutoff.Tcell, annotateTop = F, 
          genomewideline = F, suggestiveline=F)
abline(h=-log10(cutoff.Tcell), col="red")
dev.off()


# B-cell
res.dat.Bcell <- res.meta.Bcell %>% left_join(chr_bp.info, by=c("Gene" = "SetID"))
res.dat.Bcell$Chr <- as.numeric(res.dat.Bcell$Chr)
setwd("/user/xinyugao/SKATO_Project/SKATBinary_Robust_maf005/Bcell/")
png(paste0("Bcell.meta_sumz.manhattan_SKATBinary_Robust_maf005.png"), units="px", height=8000, width=15000, res=1015)
manhattan(res.dat.Bcell, chr='Chr', bp='BP', p="meta_sumz_Pvalue", snp='Gene', main="Manhattan plot of Meta Pvalues for B-cell",
          col=c("dodgerblue2","darkblue"), cex.axis=.7, las=2, ylim=c(0, 6),
          annotatePval = cutoff.Bcell, annotateTop = F, 
          genomewideline = F, suggestiveline=F)
abline(h=-log10(cutoff.Bcell), col="red")
dev.off()

# T-cell
# cohort1
res.dat.cohort1 <- Bcell.c1 %>% left_join(chr_bp.info, by = c("SetID"))
res.dat.cohort1[which.minn(res.dat.cohort1$Pvalue, 10),]
# cohort2
res.dat.cohort2 <- Bcell.c2 %>% left_join(chr_bp.info, by = c("SetID"))
res.dat.cohort2[which.minn(res.dat.cohort2$Pvalue, 10),]
# meta
res.dat.Bcell[which.minn(res.dat.Bcell$meta_sumz_Pvalue, 10),]







