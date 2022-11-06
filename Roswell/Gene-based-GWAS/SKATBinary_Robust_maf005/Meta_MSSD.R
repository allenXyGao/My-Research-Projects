mypaths <- .libPaths()
mypaths <- c(mypaths, "/projects/rpci/qzhu/xinyu/Software/R/R-4.0.2")
.libPaths(mypaths)

library(MetaSKAT)
library(dplyr)

getDiseaseCohortNames <- function(dataFile="ALLcasecontrol_cohort1") {
  split.Full_name <- unlist(strsplit(dataFile, split="_"))
  cohort.name <- split.Full_name[2]
  disease.name <- "AMLMDS"
  if (paste(unlist(strsplit(split.Full_name[1], split =""))[1:3], collapse = "") == "ALL") {
    disease.name <- "ALL"
  }
  return(list("disease_name"=disease.name, "cohort_name"=cohort.name))
}

pop <- c("ALLcasecontrol_cohort1", "ALLcasecontrol_cohort2", "AMLMDScasecontrol_cohort1", "AMLMDScasecontrol_cohort2")
get_Meta_SSDplusInfo_fromCLEANdata <- function(i) {
  names_case_cohort <- getDiseaseCohortNames(dataFile = pop[i])
  filepath.Null_Model = "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/NullModels_OrderedPheno/"
  load(file=paste0(filepath.Null_Model, "nullModel_orderedPheno_", pop[i], ".RData"))
  
  
  filepath.bin="/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/cleaned_dataset/cleaned_rare_maf005/"
  filepath.SetID="/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/SSD_results/"
  
  filepath.META_Files="/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/SKAT_META/META_SSD_RareVariants_maf005/"
  res.gen_SSD_MetaFiles <- Generate_Meta_Files(obj = res.phenotype_nullModel$null_model,
                                        
                                               File.Bed = paste0(filepath.bin, "CLEAN.rare005.", names_case_cohort$disease_name,"casecontrol_", 
                                                                 names_case_cohort$cohort_name, ".bed"),
                                               File.Bim = paste0(filepath.bin, "CLEAN.rare005.", names_case_cohort$disease_name,"casecontrol_", 
                                                                 names_case_cohort$cohort_name, ".bim"),
                                               File.SetID = paste0(filepath.SetID, "File.SetID.autosome"), 
                                               
                                               File.MSSD = paste0(filepath.META_Files, names_case_cohort$disease_name,"casecontrol_", 
                                                                  names_case_cohort$cohort_name, "/",  "File.MSSD"),
                                               File.MInfo = paste0(filepath.META_Files, names_case_cohort$disease_name,"casecontrol_", 
                                                                   names_case_cohort$cohort_name,"/",  "File.MInfo"),
                                               N.Sample = dim(res.phenotype_nullModel$phenotype_Data)[1])
  
}

get_Meta_SSDplusInfo_fromCLEANdata(1)
get_Meta_SSDplusInfo_fromCLEANdata(2)
get_Meta_SSDplusInfo_fromCLEANdata(3)
get_Meta_SSDplusInfo_fromCLEANdata(4)
