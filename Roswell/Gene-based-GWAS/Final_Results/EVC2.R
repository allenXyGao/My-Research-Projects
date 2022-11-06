### Final results for ALL 
library(SKAT)
library(MetaSKAT)

getDiseaseCohortNames <- function(dataFile="ALLcasecontrol_cohort1") {
  split.Full_name <- unlist(strsplit(dataFile, split="_"))
  cohort.name <- split.Full_name[2]
  disease.name <- "AMLMDS"
  if (paste(unlist(strsplit(split.Full_name[1], split =""))[1:3], collapse = "") == "ALL") {
    disease.name <- "ALL"
  }
  return(list("disease_name"=disease.name, "cohort_name"=cohort.name))
}


###########################################################################################################
#---------------------------------------- ALLcasecontrol_cohort1 ----------------------------------------# 
pop="ALLcasecontrol_cohort1"
names_case_cohort <- getDiseaseCohortNames(dataFile = pop)
filepath="/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/SSD_results/"
res.open_ssd <- Open_SSD(File.SSD = paste0(filepath, names_case_cohort$disease_name,"casecontrol_", names_case_cohort$cohort_name,"/", "File.SSD"),
                         File.Info = paste0(filepath, names_case_cohort$disease_name,"casecontrol_", names_case_cohort$cohort_name,"/","File.Info"))
res.open_ssd$nSets # total number of sets/Genes in this pop
res.open_ssd$SetInfo[res.open_ssd$SetInfo$SetID == "EVC2", ]
# SetIndex SetID SetSize
#   11179  EVC2      18
SetIndex <- 11179 # EVC2 index in cohort1


# pre-load phenotype data and null model
filepath.model <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/NullModels_OrderedPheno/"
load(file=paste0(filepath.model, "nullModel_orderedPheno_", pop, ".RData"))
phenotype_data <- res.phenotype_nullModel$phenotype_Data
null_model <- res.phenotype_nullModel$null_model
#pheno_with_FIDIID <- res.phenotype_nullModel$pheno_with_FIDIID

gene.name <- res.open_ssd$SetInfo[res.open_ssd$SetInfo$SetIndex == SetIndex, "SetID"]
geno.mat <- Get_Genotypes_SSD(SSD_INFO = res.open_ssd,
                              Set_Index = SetIndex,
                              is_ID = TRUE)

res.skato <- SKATBinary(geno.mat, obj = null_model, 
                        method = "optimal.adj", impute.method = "fixed")

res.skato

#---------------------------------------- ALLcasecontrol_cohort2 ----------------------------------------#
pop="ALLcasecontrol_cohort2"
names_case_cohort <- getDiseaseCohortNames(dataFile = pop)
filepath="/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/SSD_results/"
res.open_ssd <- Open_SSD(File.SSD = paste0(filepath, names_case_cohort$disease_name,"casecontrol_", names_case_cohort$cohort_name,"/", "File.SSD"),
                         File.Info = paste0(filepath, names_case_cohort$disease_name,"casecontrol_", names_case_cohort$cohort_name,"/","File.Info"))
res.open_ssd$nSets # total number of sets/Genes in this pop
res.open_ssd$SetInfo[res.open_ssd$SetInfo$SetID == "EVC2", ]
SetIndex <- 9461  # EVC2 index in cohort2
# pre-load phenotype data and null model
filepath.model <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/NullModels_OrderedPheno/"
load(file=paste0(filepath.model, "nullModel_orderedPheno_", pop, ".RData"))
phenotype_data <- res.phenotype_nullModel$phenotype_Data
null_model <- res.phenotype_nullModel$null_model
#pheno_with_FIDIID <- res.phenotype_nullModel$pheno_with_FIDIID

gene.name <- res.open_ssd$SetInfo[res.open_ssd$SetInfo$SetIndex == SetIndex, "SetID"]
geno.mat <- Get_Genotypes_SSD(SSD_INFO = res.open_ssd,
                              Set_Index = SetIndex,
                              is_ID = TRUE)

res.skato <- SKATBinary(geno.mat, obj = null_model, 
                        method.bin = "QA",
                        method = "optimal.adj", impute.method = "fixed")

res.skato


### Meta
disease <- "ALL"
filepath.cohort1 <- paste0("/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/SKAT_META/META_SSD/", 
                           disease, "casecontrol_cohort1")
filepath.cohort2 <- paste0("/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/SKAT_META/META_SSD/", 
                           disease, "casecontrol_cohort2")
File.MSSD.vec <- c(paste0(filepath.cohort1, "/File.MSSD"),
                   paste0(filepath.cohort2, "/File.MSSD"))
File.MInfo.vec <- c(paste0(filepath.cohort1, "/File.MInfo"),
                    paste0(filepath.cohort2, "/File.MInfo"))
res.open_MSSD <- Open_MSSD_File_2Read(File.MSSD.vec = File.MSSD.vec, 
                                      File.MInfo.vec = File.MInfo.vec)
temp.homo.res <- try(MetaSKAT_MSSD_OneSet(Cohort.Info = res.open_MSSD, SetID = "EVC2",
                                          method = "optimal", is.separate = FALSE))
temp.het.res <- try(MetaSKAT_MSSD_OneSet(Cohort.Info = res.open_MSSD, SetID = "EVC2",
                                         method = "optimal", is.separate = TRUE))

Info_MAF_Missing_Score_Allele12.ALL.cohort1  <- res.open_MSSD$EachInfo[[1]]$Info
Info_MAF_Missing_Score_Allele12.ALL.cohort2  <- res.open_MSSD$EachInfo[[2]]$Info
EVC2.cohort1 <- Info_MAF_Missing_Score_Allele12.ALL.cohort1[Info_MAF_Missing_Score_Allele12.ALL.cohort1$SetID=="EVC2",]
EVC2.cohort2 <- Info_MAF_Missing_Score_Allele12.ALL.cohort2[Info_MAF_Missing_Score_Allele12.ALL.cohort2$SetID=="EVC2",]
dim(EVC2.cohort1)
EVC2.cohort1
EVC2.cohort2





