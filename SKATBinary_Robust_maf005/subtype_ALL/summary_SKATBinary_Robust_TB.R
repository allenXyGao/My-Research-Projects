# summary results of batch computing
# get summarized results
library(SKAT)
library(dplyr)
library(doBy)

pop <- c("ALLcasecontrol_cohort1", "ALLcasecontrol_cohort2")
getDiseaseCohortNames <- function(dataFile="ALLcasecontrol_cohort1") {
  split.Full_name <- unlist(strsplit(dataFile, split="_"))
  cohort.name <- split.Full_name[2]
  disease.name <- "AMLMDS"
  if (paste(unlist(strsplit(split.Full_name[1], split =""))[1:3], collapse = "") == "ALL") {
    disease.name <- "ALL"
  }
  return(list("disease_name"=disease.name, "cohort_name"=cohort.name))
}


summary_SKATO_results <- function(pop, batch_size=500) {
  
  names_case_cohort <- getDiseaseCohortNames(dataFile = pop)
  filepath="/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/SSD_rare_maf005/"
  res.open_ssd <- Open_SSD(File.SSD = paste0(filepath, names_case_cohort$disease_name,"casecontrol_", names_case_cohort$cohort_name,"/", "File.SSD"),
                           File.Info = paste0(filepath, names_case_cohort$disease_name,"casecontrol_", names_case_cohort$cohort_name,"/","File.Info"))
  
  nSets.total <- res.open_ssd$nSets
  Close_SSD()
  
  start = 1; end = batch_size
  filepath.out <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/SKATO_test_results_SKATBinary_Robust_BG_maf005/"
  
  # T-cell
  out_data.Tcell <- c()
  while (start <= nSets.total) {
    end = start - 1 + batch_size
    if (end > nSets.total){
      end = nSets.total
    }
    cur_data.Tcell <- read.table(paste0(filepath.out, "Tcell/",pop, "/summary_results/", 
                                  "summary", start, "_", end, ".txt"),
                           header = TRUE)
    out_data.Tcell <- rbind(out_data.Tcell, cur_data.Tcell)
    start = end + 1
  }
  out_data.Tcell <- out_data.Tcell %>% drop_na()
  print(paste0("Summary statistics for Tcell, pop ", pop))
  print(summary(out_data.Tcell$Pvalue))
  write.table(out_data.Tcell , paste0(filepath.out, "Tcell/", pop, "/summary_results/", 
                                "allResults", ".txt"),
              row.names = FALSE, quote = FALSE)
  
  # B-cell
  out_data.Bcell <- c()
  start = 1; end = batch_size
  while (start <= nSets.total) {
    end = start - 1 + batch_size
    if (end > nSets.total){
      end = nSets.total
    }
    cur_data.Bcell <- read.table(paste0(filepath.out, "Bcell/",pop, "/summary_results/", 
                                        "summary", start, "_", end, ".txt"),
                                 header = TRUE)
    out_data.Bcell <- rbind(out_data.Bcell, cur_data.Bcell)
    start = end + 1
  }
  out_data.Bcell <- out_data.Bcell %>% drop_na()
  print(paste0("Summary statistics for B-cell, pop ", pop))
  print(summary(out_data.Bcell$Pvalue))
  write.table(out_data.Bcell , paste0(filepath.out, "Bcell/", pop, "/summary_results/", 
                                      "allResults", ".txt"),
              row.names = FALSE, quote = FALSE)

}


summary_SKATO_results(pop[1], batch_size=100)
summary_SKATO_results(pop[2], batch_size=100)


###############################################################################################



###############################################################################################
# qqplot and manhattan plot
library(qqman)
# load chr, median BP info
filepath <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/SKATO_test_results_SKATBinary_without_removing_missing/Chr_BP_Genes/"
chr_bp.info <- read.table(paste0(filepath, "chr_BP_Genes.txt"), header = TRUE)
colnames(chr_bp.info)[1] <-"SetID"
head(chr_bp.info)
# example
# head(gwasResults)
# manhattan(gwasResults)
# qq(out_data$Pvalue)

plot.fun <- function(pop) {
  names_case_cohort <- getDiseaseCohortNames(dataFile = pop)
  filepath.out <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/SKATO_test_results_SKATBinary_Robust_BG_maf005/"
  # Tcell
  d.Tcell <- read.table(paste0(filepath.out, "Tcell/", pop, "/summary_results/", 
                         "allResults", ".txt"),
                  header = TRUE)
  cutoff.Tcell = 0.05 / dim(d.Tcell)[1]
  
  # most significant genes
  print(paste0("For Tcell, pop ", pop, 
               ", the most significant gene is: ", d.Tcell[d.Tcell$Pvalue == min(d.Tcell$Pvalue), "SetID"], 
               ", with p-value: ", min(d.Tcell$Pvalue)))
  # qqplot
  setwd("/user/xinyugao/SKATO_Project/SKATBinary_Robust_maf005/Tcell/")
  png(paste0("Tcell_", pop, "SKATBinary_Robust_maf005.png")) 
  qq(d.Tcell$Pvalue, main=paste("QQ-plot of Tcell " ,  names_case_cohort$cohort_name, " with cutoff=", signif(cutoff.Tcell, 2)), cex.main=1, 
     ylim =c(0,-log10(1e-07))  )
  abline(h=-log10(cutoff.Tcell), col="blue")
  dev.off()
  # manhattan plot
  res.dat.Tcell <- d.Tcell %>% left_join(chr_bp.info, by=c("SetID"))
  res.dat.Tcell$Chr <- as.numeric(res.dat.Tcell$Chr)
  setwd("/user/xinyugao/SKATO_Project/SKATBinary_Robust_maf005/Tcell")
  png(paste0("Tcell_",pop,".manhattan_SKATBinary_Robust_maf005.png"), units="px", height=8000, width=15000, res=1015)
  manhattan(res.dat.Tcell, chr='Chr', bp='BP', p="Pvalue", snp='SetID', main=paste(" Manhattan plot of Tcell" ,  names_case_cohort$cohort_name),
            col=c("dodgerblue2","darkblue"), cex.axis=.7, las=2, ylim = c(0, 6),
            annotatePval = cutoff.Tcell, annotateTop = F, 
            genomewideline = F, suggestiveline=F)
  abline(h=-log10(cutoff.Tcell), col="red")
  dev.off()
  
  print(paste0("For Tcell, pop ", pop, ", manhattan plot is done! Cutoff value is ", cutoff.Tcell))
  
  
  # Bcell
  d.Bcell <- read.table(paste0(filepath.out, "Bcell/", pop, "/summary_results/", 
                               "allResults", ".txt"),
                        header = TRUE)
  cutoff.Bcell = 0.05 / dim(d.Bcell)[1]
  
  # most significant genes
  print(paste0("For Bcell, pop ", pop, 
               ", the most significant gene is: ", d.Bcell[d.Bcell$Pvalue == min(d.Bcell$Pvalue), "SetID"], 
               ", with p-value: ", min(d.Bcell$Pvalue)))
  # qqplot
  setwd("/user/xinyugao/SKATO_Project/SKATBinary_Robust_maf005/Bcell/")
  png(paste0("Bcell_", pop, "SKATBinary_Robust_maf005.png")) 
  qq(d.Bcell$Pvalue, main=paste("QQ-plot of Bcell " ,  names_case_cohort$cohort_name, " with cutoff=", signif(cutoff.Bcell, 2)), cex.main=1, 
     ylim =c(0,-log10(1e-07))  )
  abline(h=-log10(cutoff.Bcell), col="blue")
  dev.off()
  # manhattan plot
  res.dat.Bcell <- d.Bcell %>% left_join(chr_bp.info, by=c("SetID"))
  res.dat.Bcell$Chr <- as.numeric(res.dat.Bcell$Chr)
  setwd("/user/xinyugao/SKATO_Project/SKATBinary_Robust_maf005/Bcell")
  png(paste0("Bcell_",pop,".manhattan_SKATBinary_Robust_maf005.png"), units="px", height=8000, width=15000, res=1015)
  manhattan(res.dat.Bcell, chr='Chr', bp='BP', p="Pvalue", snp='SetID', main=paste(" Manhattan plot of Bcell" ,  names_case_cohort$cohort_name),
            col=c("dodgerblue2","darkblue"), cex.axis=.7, las=2, ylim = c(0, 6),
            annotatePval = cutoff.Bcell, annotateTop = F, 
            genomewideline = F, suggestiveline=F)
  abline(h=-log10(cutoff.Bcell), col="red")
  dev.off()
  

  print(paste0("For Bcell, pop ", pop, ", manhattan plot is done! Cutoff value is ", cutoff.Bcell))
  
}
pop <- c("ALLcasecontrol_cohort1", "ALLcasecontrol_cohort2")

# qqplot
plot.fun(pop[1])
plot.fun(pop[2])


