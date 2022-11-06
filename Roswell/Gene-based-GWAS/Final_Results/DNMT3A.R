# Final results for AMLMDS

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
#---------------------------------------- AMLMDScasecontrol_cohort1 ----------------------------------------# 
pop="AMLMDScasecontrol_cohort1"
names_case_cohort <- getDiseaseCohortNames(dataFile = pop)
filepath="/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/SSD_results/"
res.open_ssd <- Open_SSD(File.SSD = paste0(filepath, names_case_cohort$disease_name,"casecontrol_", names_case_cohort$cohort_name,"/", "File.SSD"),
                         File.Info = paste0(filepath, names_case_cohort$disease_name,"casecontrol_", names_case_cohort$cohort_name,"/","File.Info"))
res.open_ssd$nSets # total number of sets/Genes in this pop
res.open_ssd$SetInfo[res.open_ssd$SetInfo$SetID == "DNMT3A", ]
res.open_ssd$SetInfo[res.open_ssd$SetInfo$SetID == "ZG16", ]
res.open_ssd$SetInfo[res.open_ssd$SetInfo$SetID == "EXOSC5", ]
# SetIndex SetID SetSize
#   11179  EVC2      18
SetIndex_DNMT3A <- 8538 #  DNMT3A index in cohort1
SetIndex_ZG16 <- 5771 #  ZG16 index in cohort1
SetIndex_EXOSC5 <- 7997 #  EXOSC5 index in cohort1

# pre-load phenotype data and null model
filepath.model <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/NullModels_OrderedPheno/"
load(file=paste0(filepath.model, "nullModel_orderedPheno_", pop, ".RData"))
phenotype_data <- res.phenotype_nullModel$phenotype_Data
null_model <- res.phenotype_nullModel$null_model
#pheno_with_FIDIID <- res.phenotype_nullModel$pheno_with_FIDIID

# DNMT3A
gene.name <- res.open_ssd$SetInfo[res.open_ssd$SetInfo$SetIndex == SetIndex_DNMT3A, "SetID"]
geno.mat <- Get_Genotypes_SSD(SSD_INFO = res.open_ssd,
                              Set_Index = SetIndex_DNMT3A,
                              is_ID = TRUE)
res.skato <- SKATBinary(geno.mat, obj = null_model, 
                        method = "optimal.adj", impute.method = "fixed")
res.skato

# ZG16
gene.name <- res.open_ssd$SetInfo[res.open_ssd$SetInfo$SetIndex == SetIndex_ZG16, "SetID"]
geno.mat <- Get_Genotypes_SSD(SSD_INFO = res.open_ssd,
                              Set_Index = SetIndex_ZG16,
                              is_ID = TRUE)
res.skato <- SKATBinary(geno.mat, obj = null_model, 
                        method = "optimal.adj", impute.method = "fixed")
res.skato

# EXOSC5
gene.name <- res.open_ssd$SetInfo[res.open_ssd$SetInfo$SetIndex == SetIndex_EXOSC5, "SetID"]
geno.mat <- Get_Genotypes_SSD(SSD_INFO = res.open_ssd,
                              Set_Index = SetIndex_EXOSC5,
                              is_ID = TRUE)
res.skato <- SKATBinary(geno.mat, obj = null_model, 
                        method = "optimal.adj", impute.method = "fixed")
res.skato

#---------------------------------------- ALLcasecontrol_cohort2 ----------------------------------------#
pop="AMLMDScasecontrol_cohort2"
names_case_cohort <- getDiseaseCohortNames(dataFile = pop)
filepath="/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/SSD_results/"
res.open_ssd <- Open_SSD(File.SSD = paste0(filepath, names_case_cohort$disease_name,"casecontrol_", names_case_cohort$cohort_name,"/", "File.SSD"),
                         File.Info = paste0(filepath, names_case_cohort$disease_name,"casecontrol_", names_case_cohort$cohort_name,"/","File.Info"))
res.open_ssd$nSets # total number of sets/Genes in this pop
res.open_ssd$SetInfo[res.open_ssd$SetInfo$SetID == "DNMT3A", ]
res.open_ssd$SetInfo[res.open_ssd$SetInfo$SetID == "ZG16", ]
res.open_ssd$SetInfo[res.open_ssd$SetInfo$SetID == "EXOSC5", ]

SetIndex_DNMT3A <- 7661 #  DNMT3A index in cohort2
SetIndex_ZG16 <- 5175 #  ZG16 index in cohort2
SetIndex_EXOSC5 <- NA #  EXOSC5 index in cohort2

# pre-load phenotype data and null model
filepath.model <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/NullModels_OrderedPheno/"
load(file=paste0(filepath.model, "nullModel_orderedPheno_", pop, ".RData"))
phenotype_data <- res.phenotype_nullModel$phenotype_Data
null_model <- res.phenotype_nullModel$null_model
#pheno_with_FIDIID <- res.phenotype_nullModel$pheno_with_FIDIID

# DNMT3A
gene.name <- res.open_ssd$SetInfo[res.open_ssd$SetInfo$SetIndex == SetIndex_DNMT3A , "SetID"]
geno.mat <- Get_Genotypes_SSD(SSD_INFO = res.open_ssd,
                              Set_Index = SetIndex_DNMT3A ,
                              is_ID = TRUE)

res.skato <- SKATBinary(geno.mat, obj = null_model, 
                        method = "optimal.adj", impute.method = "fixed")

res.skato

# ZG16
gene.name <- res.open_ssd$SetInfo[res.open_ssd$SetInfo$SetIndex == SetIndex_ZG16 , "SetID"]
geno.mat <- Get_Genotypes_SSD(SSD_INFO = res.open_ssd,
                              Set_Index = SetIndex_ZG16 ,
                              is_ID = TRUE)

res.skato <- SKATBinary(geno.mat, obj = null_model, 
                        method = "optimal.adj", impute.method = "fixed")

res.skato




#### Meta

disease <- "AMLMDS"
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

### DNMT3A
temp.homo.res <- try(MetaSKAT_MSSD_OneSet(Cohort.Info = res.open_MSSD, SetID = "DNMT3A",
                                          method = "optimal", is.separate = FALSE))
temp.het.res <- try(MetaSKAT_MSSD_OneSet(Cohort.Info = res.open_MSSD, SetID = "DNMT3A",
                                         method = "optimal", is.separate = TRUE))

Info_MAF_Missing_Score_Allele12.cohort1  <- res.open_MSSD$EachInfo[[1]]$Info
Info_MAF_Missing_Score_Allele12.cohort2  <- res.open_MSSD$EachInfo[[2]]$Info
DNMT3A.cohort1 <- Info_MAF_Missing_Score_Allele12.cohort1[Info_MAF_Missing_Score_Allele12.cohort1$SetID=="DNMT3A",]
DNMT3A.cohort2 <- Info_MAF_Missing_Score_Allele12.cohort2[Info_MAF_Missing_Score_Allele12.cohort2$SetID=="DNMT3A",]
DNMT3A.cohort1
DNMT3A.cohort2

### ZG16
temp.homo.res <- try(MetaSKAT_MSSD_OneSet(Cohort.Info = res.open_MSSD, SetID = "ZG16",
                                          method = "optimal", is.separate = FALSE))
temp.het.res <- try(MetaSKAT_MSSD_OneSet(Cohort.Info = res.open_MSSD, SetID = "ZG16",
                                         method = "optimal", is.separate = TRUE))

Info_MAF_Missing_Score_Allele12.cohort1  <- res.open_MSSD$EachInfo[[1]]$Info
Info_MAF_Missing_Score_Allele12.cohort2  <- res.open_MSSD$EachInfo[[2]]$Info
DNMT3A.cohort1 <- Info_MAF_Missing_Score_Allele12.cohort1[Info_MAF_Missing_Score_Allele12.cohort1$SetID=="ZG16",]
DNMT3A.cohort2 <- Info_MAF_Missing_Score_Allele12.cohort2[Info_MAF_Missing_Score_Allele12.cohort2$SetID=="ZG16",]
DNMT3A.cohort1
DNMT3A.cohort2


### EXOSC5
temp.homo.res <- try(MetaSKAT_MSSD_OneSet(Cohort.Info = res.open_MSSD, SetID = "EXOSC5",
                                          method = "optimal", is.separate = FALSE))
temp.het.res <- try(MetaSKAT_MSSD_OneSet(Cohort.Info = res.open_MSSD, SetID = "EXOSC5",
                                         method = "optimal", is.separate = TRUE))

Info_MAF_Missing_Score_Allele12.cohort1  <- res.open_MSSD$EachInfo[[1]]$Info
Info_MAF_Missing_Score_Allele12.cohort2  <- res.open_MSSD$EachInfo[[2]]$Info
DNMT3A.cohort1 <- Info_MAF_Missing_Score_Allele12.cohort1[Info_MAF_Missing_Score_Allele12.cohort1$SetID=="EXOSC5",]
DNMT3A.cohort2 <- Info_MAF_Missing_Score_Allele12.cohort2[Info_MAF_Missing_Score_Allele12.cohort2$SetID=="EXOSC5",]
DNMT3A.cohort1
DNMT3A.cohort2
