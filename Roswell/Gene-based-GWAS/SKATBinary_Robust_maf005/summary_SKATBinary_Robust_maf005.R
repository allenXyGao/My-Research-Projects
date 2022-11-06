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

  out_data <- c()
  while (start <= nSets.total) {
    end = start - 1 + batch_size
    if (end > nSets.total){
      end = nSets.total
    }
    cur_data <- read.table(paste0(filepath.out, pop, "/summary_results/", 
                                  "summary", start, "_", end, ".txt"),
                           header = TRUE)
    out_data <- rbind(out_data, cur_data)
    start = end + 1
  }
  print(paste0("For pop ", pop, ", ", nSets.total - dim(out_data)[1],
               " genes with less than 2 variants have been removed"))
  print(paste0("Summary statistics for pop ", pop))
  print(summary(out_data$Pvalue))
  
  write.table(out_data , paste0(filepath.out, pop, "/summary_results/", 
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
  d <- read.table(paste0(filepath.out, pop, "/summary_results/", 
                         "allResults", ".txt"),
                  header = TRUE)
  
  
  cutoff = 0.05 / dim(d)[1]
  
  
  # most significant genes
  print(paste0("For pop ", pop, 
               ", the most significant gene is: ", d[d$Pvalue == min(d$Pvalue), "SetID"], 
               ", with p-value: ", min(d$Pvalue)))
  
  # qqplot
  setwd("/user/xinyugao/SKATO_Project/SKATBinary_Robust_maf005/")
  png(paste0(pop, "SKATBinary_Robust_maf005.png")) 
  qq(d$Pvalue, main=paste("QQ-plot of " ,  pop, " with cutoff=", signif(cutoff, 2)), cex.main=1, 
     ylim =c(0,-log10(1e-07))  )
  abline(h=-log10(cutoff), col="blue")
  dev.off()
  
  

  # manhattan plot
  res.dat <- d %>% left_join(chr_bp.info, by=c("SetID"))
  res.dat$Chr <- as.numeric(res.dat$Chr)
  setwd("/user/xinyugao/SKATO_Project/SKATBinary_Robust_maf005/")
  png(paste0(pop,".manhattan_SKATBinary_Robust_maf005.png"), units="px", height=8000, width=15000, res=1015)
  manhattan(res.dat, chr='Chr', bp='BP', p="Pvalue", snp='SetID', main=paste(pop),
            col=c("dodgerblue2","darkblue"), cex.axis=.7, las=2, 
            annotatePval = cutoff, annotateTop = F, 
            genomewideline = F, suggestiveline=F)
  abline(h=-log10(cutoff), col="red")
  dev.off()
  print(paste0("For pop ", pop, ", manhattan plot is done! Cutoff value is ", cutoff))
  
}
pop <- c("ALLcasecontrol_cohort1", "ALLcasecontrol_cohort2")

# qqplot
plot.fun(pop[1])
plot.fun(pop[2])


