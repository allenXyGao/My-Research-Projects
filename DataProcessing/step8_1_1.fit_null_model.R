# this R script is used to fit null models for each pop
library(SKAT)

pop <- c("ALLcasecontrol_cohort1", "ALLcasecontrol_cohort2",
         "AMLMDScasecontrol_cohort1", "AMLMDScasecontrol_cohort2")
getDiseaseCohortNames <- function(dataFile="ALLcasecontrol_cohort1") {
  split.Full_name <- unlist(strsplit(dataFile, split="_"))
  cohort.name <- split.Full_name[2]
  disease.name <- "AMLMDS"
  if (paste(unlist(strsplit(split.Full_name[1], split =""))[1:3], collapse = "") == "ALL") {
    disease.name <- "ALL"
  }
  return(list("disease_name"=disease.name, "cohort_name"=cohort.name))
}

########################################
### Distribution of "gender" in each pop:
# ALLcasecontrol_cohort1: 
#       Female   Male
#          803   1622
# ALLcasecontrol_cohort2
#       Female   Male
#          139  416
# AMLMDScasecontrol_cohort1
#       Female   Male
#         1331   2179
# AMLMDScasecontrol_cohort2

#       Female   Male
#          332   642
###################################################################
#### Distribution of the binary trait "disease" or "non-disease"
# ALLcasecontrol_cohort1: 
#       case         control
#     442 (18.2%)    1983 (81.8%)
# ALLcasecontrol_cohort2
#       case         control
#     41 (7.4%)      514 (92.6%)
# AMLMDScasecontrol_cohort1
#       case         control
#    1525 (43.4%)    1985 (56.6%)
# AMLMDScasecontrol_cohort2
#       case         control
#     460 (47.2%)    514 (52.8%)
###################################################################

get_nullModel_phenotypeData <- function(i) {
  ################################################################################
  # This function aims to create phenotype datasets with the same order
  # of FAM file, so that the sample in both genotype matrix and null_model
  # can be matched!!
  ################################################################################
  
  names_case_cohort <- getDiseaseCohortNames(pop[i])
  filepath.pheno <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/disorder_Penotype_Data/phenotypeData_without_Step7_2/"
  pheno.data <- read.table(paste0(filepath.pheno, "phenotype.data.", names_case_cohort$disease_name,
                                  ".casecontrol.", names_case_cohort$cohort_name),
                           header = T)
  
  pheno.data[pheno.data$sex == "Male", "sex"] <- 1
  pheno.data[pheno.data$sex == "Female", "sex"] <- 0
  pheno.data$sex <- as.factor(pheno.data$sex)
  pheno.data[pheno.data$phenotype == 1, "phenotype"] <- 0 # control
  pheno.data[pheno.data$phenotype == 2, "phenotype"] <- 1 # case
  
  ################################################################################
  filepath.bin="/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/cleaned_dataset/cleaned_rare_variants/"
  File.Fam = paste0(filepath.bin, "CLEAN.rare.", names_case_cohort$disease_name,"casecontrol_", names_case_cohort$cohort_name, ".fam")
  FAM <- Read_Plink_FAM(File.Fam)
  # recover order
  pheno.data <- pheno.data[match(FAM$IID, pheno.data$IID), ]
  if (sum(pheno.data$FID != FAM$FID) > 0 ) {
    print(paste0("warning: For pop ", pop[i], "phenotype data has mismatch order!!!"))
  }
  ################################################################################

  data_pheno_covariates <- pheno.data[, -c(1,2)]
  
  null_model <- SKAT_Null_Model(phenotype ~ .,
                                data = data_pheno_covariates , out_type = "D")
  filepath.model <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/NullModels_OrderedPheno/"
  
  # store odered phenotype data and null_model
  res.phenotype_nullModel <- list()
  res.phenotype_nullModel$phenotype_Data <- data_pheno_covariates
  res.phenotype_nullModel$pheno_with_FIDIID <- pheno.data
  res.phenotype_nullModel$null_model <- null_model
  save(res.phenotype_nullModel, file=paste0(filepath.model,"nullModel_orderedPheno_", pop[i],".RData"))
  print(paste0("For pop ", pop[i], ", Done!"))
}

get_nullModel_phenotypeData(1)
get_nullModel_phenotypeData(2)
get_nullModel_phenotypeData(3)
get_nullModel_phenotypeData(4)


# e.g.
rm(list=ls())
i=1
pop <- c("ALLcasecontrol_cohort1", "ALLcasecontrol_cohort2",
         "AMLMDScasecontrol_cohort1", "AMLMDScasecontrol_cohort2")
filepath.model <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/NullModels_OrderedPheno/"
load(file=paste0(filepath.model, "nullModel_orderedPheno_", pop[i], ".RData"))

table(res.phenotype_nullModel$phenotype_Data$phenotype)
prop.table(table(res.phenotype_nullModel$phenotype_Data$phenotype))
