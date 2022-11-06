# auto_gen_Phenotype_DataFiles

rm(list = ls())
library(dplyr)
library(tidyr)

pop <- c("ALLcasecontrol_cohort1", "ALLcasecontrol_cohort2",
         "AMLMDScasecontrol_cohort1", "AMLMDScasecontrol_cohort2")
disease <- c("ALL", "AML.MDS")

########################### load and clean meta_phenotype data ########################
filepath <- "/projects/rpci/qzhu/lsuchest/lsuchest/DBMT_PhenoData/"
#meta_phenotype <- read.table(file=paste(filepath, "DBMT_PhenoData_EA_long_allVar_20190223.txt", sep=""), fill = T)
#meta_phenotype <- read.table(file=paste(filepath, "DBMT_PhenoData_EA_long_allVar_20190223.txt", sep=""), sep="\t", header = T)
# /projects/rpci/qzhu/lsuchest/lsuchest/DBMT_PhenoData/DBMT_PhenoData_long_allVar_20180424.txt
meta_phenotype <- read.table(file=paste(filepath, "DBMT_PhenoData_long_allVar_20180424.txt", sep=""), sep="\t", header = T)
colnames(meta_phenotype)
meta_phenotype.less <- meta_phenotype[, c("Exome.FID", "Exome.IID", "Exome.recip.pop", "Exome.donor.pop", 
                                          "age", "dnrage", "sex", "dnrsex")]
meta_phenotype.less$Exome.FID <- as.character(meta_phenotype.less$Exome.FID)
meta_phenotype.less$Exome.IID <- as.character(meta_phenotype.less$Exome.IID)

colnames(meta_phenotype.less)[colnames(meta_phenotype.less) == "Exome.FID"] <- "FID"
colnames(meta_phenotype.less)[colnames(meta_phenotype.less) == "Exome.IID"] <- "IID"
head(meta_phenotype.less)
dim(meta_phenotype.less)
colSums((!is.na(meta_phenotype.less)))
#######################################################################################

gen_data <- function(i, merge_method="inner") {
  #####################################################################
  # i: index of 4 datafiles, where 1:ALLcasecontrol_cohort1, 2:ALLcasecontrol_cohort2
  #                                3:AMLMDScasecontrol_cohort1, 4:AMLMDScasecontrol_cohort2
  # merge_method: merge meta_phenotype dataset and PCs+phenotype dataset
  # (since not all samples can be matched)
  #####################################################################
  
  #-----------------------------step 1: clean PC2 file and case/control file-----------------------------
  # PC
  filepath <- paste("/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/EigenStrat_Results/",
                    pop[i],"/", sep="")
  PCs <- read.table(file=paste(filepath, "LDprune_autosome-pca2.eth.out", sep=""), header = T)
  PC.names <- c()
  for (k in 1:20) {
    PC.names <- c(PC.names, paste0("PC", k))
  }
  colnames(PCs)[3:dim(PCs)[2]] <- PC.names
  # extract significant PCs for this specific pop
  twout <- read.table(file=paste(filepath, "LDprune_autosome-pca/LDprune_autosome-pca2.twout", sep=""))
  sig.n <- sum(twout[1:20,5]<=0.05)  ## use top 20 to calculate, sig.n<20
  print(paste("The number of significant PCs=", sig.n))
  PCs.sig.cut <- PCs[, c(3:(3+sig.n-1))]
  #head(PCs.sig.cut)
  # identify FID/IID
  out <- strsplit(as.character(PCs[, 1]),':')
  FID_IID.remaining <- do.call(rbind, out)
  PCs.sig.FID.IID <- cbind(FID_IID.remaining, PCs.sig.cut)
  colnames(PCs.sig.FID.IID)[1:2] <- c("FID", "short_IID")
  #print(head(PCs.sig.FID.IID))
  
  ## case/control
  # case=2
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
  case.data$phenotype <- 2
  # control=1
  control.filepath <- paste("/projects/rpci/qzhu/qzhu/mywork/LaraBMTexome/input/BMT.snpQC.final.donor.EA.",
                            cohort.name, ".fam", sep="")
  control.data <- read.table(file=control.filepath)[, 1:2]
  control.data$phenotype <- 1
  casecontrol.data <- rbind(case.data, control.data)
  colnames(casecontrol.data)[1:2] <- c("FID", "IID")
 
  
  #-----------------------------step 2: merge datafiles-----------------------------
  # merge PCs.sig.FID.IID and ALL.casecontrol.cohort1 -> res1
  out2 <- strsplit(casecontrol.data$IID, "([A-Z]-)", perl = TRUE)
  casecontrol.data$short_IID <- do.call(rbind, out2)[,2]
  casecontrol.data$FID <- as.character(casecontrol.data$FID)
  casecontrol.data$short_IID <- as.character(casecontrol.data$short_IID)
  PCs.sig.FID.IID$FID <- as.character(PCs.sig.FID.IID$FID)
  PCs.sig.FID.IID$short_IID <- as.character(PCs.sig.FID.IID$short_IID)
  res1 <- casecontrol.data %>% inner_join(PCs.sig.FID.IID, by=c("FID", "short_IID"))
  #print(head(res1))
  
  # merge res1 and meta_phenotype.less
  if (merge_method == "inner") {
    res2 <- res1 %>% inner_join(meta_phenotype.less, by=c("FID", "IID"))
    #res2 <- res1 %>% inner_join(meta_phenotype.less, by="IID")
  }
  if (merge_method == "left") {
    res2 <- res1 %>% left_join(meta_phenotype.less, by=c("FID", "IID"))
  }
  
  #-----------------------------step 3: select columns-----------------------------
  # cases=2, control=1
  res2$data.sex <- NA
  res2$data.sex[res2$phenotype == 1] <-  res2[res2$phenotype == 1, "dnrsex"]
  res2$data.sex[res2$phenotype == 2] <-  res2[res2$phenotype == 2, "sex"]
  res2$data.age <- NA
  res2$data.age[res2$phenotype == 1] <- res2[res2$phenotype == 1, "dnrage"]
  res2$data.age[res2$phenotype == 2] <- res2[res2$phenotype == 2, "age"]
  
  res.phenotype <- res2[, c("FID", "IID", "phenotype", "data.age", "data.sex", PC.names[1:sig.n])]
  colnames(res.phenotype)[colnames(res.phenotype) == "data.age"] = "age"
  colnames(res.phenotype)[colnames(res.phenotype) == "data.sex"] = "sex"
  
  
  return(res.phenotype)
  
  
}


#-------------------------------------Results-------------------------------------------
# 1:ALLcasecontrol_cohort1,    2:ALLcasecontrol_cohort2
# 3:AMLMDScasecontrol_cohort1, 4:AMLMDScasecontrol_cohort2

## ALL case control cohort1 
# ALL casecontrol cohort1 dataset using inner join
phenotype.data.ALL.casecontrol.cohort1.inner <- gen_data(1, merge_method = "inner")
dim(phenotype.data.ALL.casecontrol.cohort1.inner)
head(phenotype.data.ALL.casecontrol.cohort1.inner)
# ALL casecontrol cohort1 dataset using left join
phenotype.data.ALL.casecontrol.cohort1.left <- gen_data(1, merge_method = "left")
head(phenotype.data.ALL.casecontrol.cohort1.left)
dim(phenotype.data.ALL.casecontrol.cohort1.left)

## ALL case control cohort2
# ALL casecontrol cohort2 dataset using inner join
phenotype.data.ALL.casecontrol.cohort2.inner <- gen_data(2, merge_method = "inner")
dim(phenotype.data.ALL.casecontrol.cohort2.inner)
head(phenotype.data.ALL.casecontrol.cohort2.inner)
# ALL casecontrol cohort1 dataset using left join
phenotype.data.ALL.casecontrol.cohort2.left <- gen_data(2, merge_method = "left")
head(phenotype.data.ALL.casecontrol.cohort2.left)
dim(phenotype.data.ALL.casecontrol.cohort2.left)

## AMS.MDS case control cohort1
# AMS.MDS casecontrol cohort1 dataset using inner join
phenotype.data.AMS.MDS.casecontrol.cohort1.inner <- gen_data(3, merge_method = "inner")
head(phenotype.data.AMS.MDS.casecontrol.cohort1.inner)
dim(phenotype.data.AMS.MDS.casecontrol.cohort1.inner)
# AMS.MDS casecontrol cohort1 dataset using left join
phenotype.data.AMS.MDS.casecontrol.cohort1.left <- gen_data(3, merge_method = "left")
head(phenotype.data.AMS.MDS.casecontrol.cohort1.left)
dim(phenotype.data.AMS.MDS.casecontrol.cohort1.left)


## AMS.MDS case control cohort2
# AMS.MDS casecontrol cohort2 dataset using inner join
phenotype.data.AMS.MDS.casecontrol.cohort2.inner <- gen_data(4, merge_method = "inner")
head(phenotype.data.AMS.MDS.casecontrol.cohort2.inner)
dim(phenotype.data.AMS.MDS.casecontrol.cohort2.inner)
# AMS.MDS casecontrol cohort2 dataset using left join
phenotype.data.AMS.MDS.casecontrol.cohort2.left <- gen_data(4, merge_method = "left")
head(phenotype.data.AMS.MDS.casecontrol.cohort2.left)
dim(phenotype.data.AMS.MDS.casecontrol.cohort2.left)




#-----------------------------------Store Data-------------------------------------------
# filepath="/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/Penotype_Data/phenotypeData_without_Step7_2"
# write.table(genotype.data.ALL.casecontrol.cohort1.inner, 
#             file = paste0(filepath, "/", "genotype.data.ALL.casecontrol.cohort1"), 
#             sep="\t")

store_data <- function(i) {
  filepath="/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/Penotype_Data/phenotypeData_without_Step7_2"
  temp_data <- gen_data(i, merge_method = "inner")
  split.Full_name <- unlist(strsplit(pop[i], split="_"))
  cohort.name <- split.Full_name[2]
  disease.name <- "AML.MDS"
  if (paste(unlist(strsplit(split.Full_name[1], split =""))[1:3], collapse = "") == "ALL") {
    disease.name <- "ALL"
  }
  # store data into txt files
  write.table(temp_data, 
    file = paste0(filepath, "/", "phenotype.data.", disease.name, ".casecontrol.", cohort.name),
    sep="\t", row.names = FALSE)
  print("Done")
}

sapply(1:4, store_data)



#-----------------------------------Double Check-----------------------------------------
## check by /fam file
# e.g. /projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/ALLcasecontrol_cohort2.fam
pop <- c("ALLcasecontrol_cohort1", "ALLcasecontrol_cohort2",
         "AMLMDScasecontrol_cohort1", "AMLMDScasecontrol_cohort2")
disease <- c("ALL", "AML.MDS")

checkSexWithFam <- function(i) {
  fam_filepath="/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/"
  pheno_filepath="/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/Penotype_Data/phenotypeData_without_Step7_2/"
  split.Full_name <- unlist(strsplit(pop[i], split="_"))
  cohort.name <- split.Full_name[2]
  disease.name.pheno <- "AML.MDS"
  disease.name.fam <- "AMLMDS"
  if (paste(unlist(strsplit(split.Full_name[1], split =""))[1:3], collapse = "") == "ALL") {
    disease.name.pheno <- "ALL"
    disease.name.fam <- "ALL"
  }
  pheno_temp <- read.table(file=paste0(pheno_filepath, "phenotype.data.", disease.name.pheno, ".casecontrol.", cohort.name), header = T)
  fam_temp <- read.table(file=paste0(fam_filepath, disease.name.fam, "casecontrol_",cohort.name, ".fam"), header = F)
  colnames(fam_temp) <- c("FID", "IID", "Father", "Mother", "Gender", "Pheno")
  fam_temp[fam_temp$Gender == 1, "Gender"] <- "Male"
  fam_temp[fam_temp$Gender == 2, "Gender"] <- "Female"
  res.compare <- pheno_temp %>% inner_join(fam_temp, by=c("FID", "IID"))
  
  if (sum(res.compare$sex != res.compare$Gender) == 0) {
    print(paste("For",pop[i], ", Sex Check passed!"))
  }
  else {
    print(paste("For",pop[i], ", Sex Check Failed!"))
  }
  
}

# check pass/fail
checkSexWithFam(1)
checkSexWithFam(2)
checkSexWithFam(3)
checkSexWithFam(4)

