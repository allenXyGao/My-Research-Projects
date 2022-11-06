# Step 1: Removed samples for each pop
library(dplyr)
library(tidyr)

pop <- c("ALLcasecontrol_cohort1", "ALLcasecontrol_cohort2",
        "AMLMDScasecontrol_cohort1", "AMLMDScasecontrol_cohort2")

get_removed_samples_from_pca2 <- function(i) {
  filepath="/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/EigenStrat_Results/"
  print(pop[i])
  temp_removed <- read.table(paste0(filepath, pop[i], "/LDprune_autosome-pca/","removedSamples.", pop[i]))
  colnames(temp_removed) <-  c("FID", "short_IID")
  # recover its IID
  split.Full_name <- unlist(strsplit(pop[i], split="_"))
  cohort.name <- split.Full_name[2]
  disease.name <- "AML.MDS"
  if (paste(unlist(strsplit(split.Full_name[1], split =""))[1:3], collapse = "") == "ALL") {
    disease.name <- "ALL"
  }
  case.filepath <- paste("/projects/rpci/qzhu/eschille/transplant.GRANT/data/Recipients/EA",
                         cohort.name,  "/", 
                         disease.name, ".samples.txt", sep="")
  
  case.data <- read.table(file=case.filepath)
  control.filepath <- paste("/projects/rpci/qzhu/qzhu/mywork/LaraBMTexome/input/BMT.snpQC.final.donor.EA.",
                            cohort.name, ".fam", sep="")
  control.data <- read.table(file=control.filepath)[, 1:2]
  casecontrol.data <- rbind(case.data, control.data)
  colnames(casecontrol.data)[1:2] <- c("FID", "IID")
  
  out <- strsplit(casecontrol.data$IID, "([A-Z]-)", perl = TRUE)
  casecontrol.data$short_IID <- do.call(rbind, out)[,2]
  casecontrol.data$FID <- as.character(casecontrol.data$FID)
  casecontrol.data$short_IID <- as.character(casecontrol.data$short_IID)
  temp_removed$FID <- as.character(temp_removed$FID)
  temp_removed$short_IID <- as.character(temp_removed$short_IID)
  res <- casecontrol.data %>% inner_join(temp_removed, by=c("FID", "short_IID"))
  res <- res[, c(1,2)]
  filepath.out = "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/removed_Samples_SNPs"
  disease.name <- "AMLMDS"
  if (paste(unlist(strsplit(split.Full_name[1], split =""))[1:3], collapse = "") == "ALL") {
    disease.name <- "ALL"
  }
  write.table(res, 
              file = paste0(filepath.out, "/", "removed.samples_", 
                            disease.name, "casecontrol_", cohort.name),
              sep=" ", row.names = FALSE,  quote = FALSE)
  
  return("Recover IID success")
  
}

get_removed_samples_from_pca2(1)
get_removed_samples_from_pca2(2)
get_removed_samples_from_pca2(3)
get_removed_samples_from_pca2(4)




#----------------------------------------------------------------------------------#
# Step 2: combine all excluded SNP lists (e.g. 0/0) to All_excludeSNPs.txt
# please check the original .bim files you created to see whether there are any
# SNPs with both alleles 0 (the 5th and 6th column)
# in .bim files, column 1: chromsome code; column 2. variant identifier;
#                column 5: Allele 1 (usually minor); column 6. Allel 2 (usually major)

pop <- c("ALLcasecontrol_cohort1", "ALLcasecontrol_cohort2",
         "AMLMDScasecontrol_cohort1", "AMLMDScasecontrol_cohort2")

getExcludedSNPs <- function(i) {
  filepath <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/"
  temp_bim <- read.table(paste0(filepath, pop[i], ".bim"))
  # print(pop[i])
  # print(head(temp_bim))
  # print(dim(temp_bim))
  excludedSNPs <- temp_bim[(temp_bim$V5 == 0 && temp_bim$V6 == 0), "V2"]
  return(excludedSNPs)
}
getExcludedSNPs(1)
getExcludedSNPs(2)
getExcludedSNPs(3)
getExcludedSNPs(4)
# No SNPs removed

# check 
checkFreqAllele <- function(i) {
  filepath <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/"
  temp_bim <- read.table(paste0(filepath, pop[i], ".bim"))
  print(paste("#----------", pop[i], "------------#"))
  print(paste0("numbers of SNPs with column 5 == 0, ", dim(temp_bim[temp_bim$V5 == 0,])[1]))
  print(paste0("numbers of SNPs with column 6 == 0, ", dim(temp_bim[temp_bim$V6 == 0,])[1]))
  print("Freq of allele in column 5")
  print(table(temp_bim$V5))
  print("----------------------------------------------")
  print("Freq of allele in column 6")
  print(table(temp_bim$V6))
  print("#------------------Done-----------------#")
}
checkFreqAllele(1)
checkFreqAllele(2)
checkFreqAllele(3)
checkFreqAllele(4)
