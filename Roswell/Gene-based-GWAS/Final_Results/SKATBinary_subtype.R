# T-ALL & B-ALL SKATBinary_robust
rm(list = ls())
library(SKAT)

# cohort1
filepath <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/subtypes_Tcell_Bcell/"
load(file=paste0(filepath, "ALL.EVC2.cohort1.subtype.info.RData"))
case.Tcell.IID.cohort1 <- c(ALL.EVC2.cohort1.subtype.info$carriers_TCell_IID, ALL.EVC2.cohort1.subtype.info$noncarriers_TCell_IID)
case.Bcell.IID.cohort1 <- c(ALL.EVC2.cohort1.subtype.info$carriers_BCell_IID, ALL.EVC2.cohort1.subtype.info$noncarriers_BCell_ID)

#---------cohort1-------------#
# subtype
# carriers T-cell B-cell  (unknow)  Control
# YES      5       21         6        47
# NO      58      289        63      1936
# total   63      310        69      1983

# T-cell
# case control
# 63    1983
# 3%    97%
# length(case.Tcell.IID.cohort1)
# B-cell
# case  control
# 310    1983
# 13.5%   86.5%
# length(case.Bcell.IID.cohort1)



# cohort2
filepath <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/subtypes_Tcell_Bcell/"
load(file=paste0(filepath, "ALL.EVC2.cohort2.subtype.info.RData"))
case.Tcell.IID.cohort2 <- c(ALL.EVC2.cohort2.subtype.info$carriers_TCell_IID, ALL.EVC2.cohort2.subtype.info$noncarriers_TCell_IID)
case.Bcell.IID.cohort2 <- c(ALL.EVC2.cohort2.subtype.info$carriers_BCell_IID, ALL.EVC2.cohort2.subtype.info$noncarriers_BCell_ID)

#---------cohort2-------------#
# subtype
# carriers T-cell B-cell  (unknow)  Control
# YES      0        3         0       17
# NO       4       32         2       497
# total    4       35         2       514

# T-cell
# case   control
#  4       514
# 0.77%   99.23%
# B-cell
# case   control
# 35      514
# 6.38%   93.62%




getDiseaseCohortNames <- function(dataFile="ALLcasecontrol_cohort1") {
  split.Full_name <- unlist(strsplit(dataFile, split="_"))
  cohort.name <- split.Full_name[2]
  disease.name <- "AMLMDS"
  if (paste(unlist(strsplit(split.Full_name[1], split =""))[1:3], collapse = "") == "ALL") {
    disease.name <- "ALL"
  }
  return(list("disease_name"=disease.name, "cohort_name"=cohort.name))
}


subtype.SKATBinary.robust <- function(pop= "ALLcasecontrol_cohort1", gene.name="EVC2", impute.method="bestguess") {
  # -----------------------------load some necessary info -------------------------------------#
  names_case_cohort <- getDiseaseCohortNames(dataFile = pop)
  if (names_case_cohort$cohort_name == "cohort1") {
    case.Tcell.IID <- case.Tcell.IID.cohort1; case.Bcell.IID <- case.Bcell.IID.cohort1
  }
  else if (names_case_cohort$cohort_name == "cohort2") {
    case.Tcell.IID <- case.Tcell.IID.cohort2; case.Bcell.IID <- case.Bcell.IID.cohort2
  }
  
  filepath="/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/SSD_results/"
  res.open_ssd <- Open_SSD(File.SSD = paste0(filepath, names_case_cohort$disease_name,"casecontrol_", names_case_cohort$cohort_name,"/", "File.SSD"),
                           File.Info = paste0(filepath, names_case_cohort$disease_name,"casecontrol_", names_case_cohort$cohort_name,"/","File.Info"))
  
  filepath.model <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/NullModels_OrderedPheno/"
  load(file=paste0(filepath.model, "nullModel_orderedPheno_", pop, ".RData"))
  phenotype_data <- res.phenotype_nullModel$phenotype_Data
  pheno_with_FIDIID <- res.phenotype_nullModel$pheno_with_FIDIID
  index <- res.open_ssd$SetInfo[res.open_ssd$SetInfo$SetID == gene.name, "SetIndex"]
  geno.mat <- Get_Genotypes_SSD(SSD_INFO = res.open_ssd,
                                Set_Index = index, is_ID = TRUE)
  # --------------------------------------------------------------------------------------------#
  
  #--------------------------------T-cell SKATBinary_robust------------------------------------------
  control.IID <- pheno_with_FIDIID$IID[pheno_with_FIDIID$phenotype == 0]
  Tcell.casecontrol.IID <- c(case.Tcell.IID, control.IID)
  # phenotype data
  casecontrol.Tcell.dat <- phenotype_data[!is.na(match(pheno_with_FIDIID$IID, Tcell.casecontrol.IID)), ]
  Tcell.sampleList <- pheno_with_FIDIID[!is.na(match(pheno_with_FIDIID$IID, Tcell.casecontrol.IID)), c("FID", "IID")]
  # genotype data
  geno.mat.Tcell <- geno.mat[!is.na(match(pheno_with_FIDIID$IID, Tcell.casecontrol.IID)), ]
  # refit-null model
  model.Tcell <- SKAT_Null_Model(phenotype ~ .,
                           data = casecontrol.Tcell.dat , out_type = "D")
  # flag <- apply(geno.mat.Tcell , 1, function(x) {any(is.na(x) | x == 9)})
  res.skato.Tcell <- SKATBinary_Robust(geno.mat.Tcell, obj = model.Tcell, method = "SKATO")
  #--------------------------------------------------------------------------------------------------
  
  #--------------------------------B-cell SKATBinary_robust------------------------------------------
  Bcell.casecontrol.IID <- c(case.Bcell.IID, control.IID)
  # phenotype data
  casecontrol.Bcell.dat <- phenotype_data[!is.na(match(pheno_with_FIDIID$IID, Bcell.casecontrol.IID)), ]
  Bcell.sampleList <- pheno_with_FIDIID[!is.na(match(pheno_with_FIDIID$IID, Bcell.casecontrol.IID)), c("FID", "IID")]
  # genotype data
  geno.mat.Bcell <- geno.mat[!is.na(match(pheno_with_FIDIID$IID, Bcell.casecontrol.IID)), ]
  # refit-null model
  model.Bcell <- SKAT_Null_Model(phenotype ~ .,
                           data = casecontrol.Bcell.dat , out_type = "D")
  # flag <- apply(geno.mat.Tcell , 1, function(x) {any(is.na(x) | x == 9)})
  res.skato.Bcell <- SKATBinary_Robust(geno.mat.Bcell, obj = model.Bcell, method = "SKATO", impute.method = impute.method)
  #--------------------------------------------------------------------------------------------------
  
  return(list("Tcell.res"=res.skato.Tcell, "Bcell.res"=res.skato.Bcell, 
              "nullModel.Tcell"=model.Tcell, "nullModel.Bcell"=model.Bcell,
              "Tcell.casecontrol.sampleList"=Tcell.sampleList, "Bcell.casecontrol.sampleList"=Bcell.sampleList))
  
}


subtype.res.cohort1 <- subtype.SKATBinary.robust(pop= "ALLcasecontrol_cohort1", gene.name="EVC2")
subtype.res.cohort2 <- subtype.SKATBinary.robust(pop= "ALLcasecontrol_cohort2", gene.name="EVC2")

# store results
# just store sample List
filepath.out <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/subtypes_Tcell_Bcell/Tcell_Meta/"
write.table(subtype.res.cohort1$Tcell.casecontrol.sampleList, file = paste0(filepath.out, "ALLEVC2.cohort1.Tcell.sampleList.txt"),
            col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(subtype.res.cohort2$Tcell.casecontrol.sampleList, file = paste0(filepath.out, "ALLEVC2.cohort2.Tcell.sampleList.txt"),
            col.names = FALSE, row.names = FALSE, quote = FALSE)

filepath.out <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/subtypes_Tcell_Bcell/Bcell_Meta/"
write.table(subtype.res.cohort1$Bcell.casecontrol.sampleList, file = paste0(filepath.out, "ALLEVC2.cohort1.Bcell.sampleList.txt"),
            col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(subtype.res.cohort2$Bcell.casecontrol.sampleList, file = paste0(filepath.out, "ALLEVC2.cohort2.Bcell.sampleList.txt"),
            col.names = FALSE, row.names = FALSE, quote = FALSE)
# store all results
filepath.out <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/subtypes_Tcell_Bcell/"
save(subtype.res.cohort1, file=paste0(filepath.out, "ALL.EVC2.cohort1_nullModel_testRes_sampleList.RData"))
save(subtype.res.cohort2, file=paste0(filepath.out, "ALL.EVC2.cohort2_nullModel_testRes_sampleList.RData"))

  
  
# compare with all cases versus control using SKATBinary_robust

ALLcase.SKATBinary_robust <- function(pop= "ALLcasecontrol_cohort1", gene.name="EVC2", impute.method="bestguess") {
  names_case_cohort <- getDiseaseCohortNames(dataFile = pop)
  filepath="/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/SSD_results/"
  res.open_ssd <- Open_SSD(File.SSD = paste0(filepath, names_case_cohort$disease_name,"casecontrol_", names_case_cohort$cohort_name,"/", "File.SSD"),
                           File.Info = paste0(filepath, names_case_cohort$disease_name,"casecontrol_", names_case_cohort$cohort_name,"/","File.Info"))
  filepath.model <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/NullModels_OrderedPheno/"
  load(file=paste0(filepath.model, "nullModel_orderedPheno_", pop, ".RData"))
  model <- res.phenotype_nullModel$null_model
  index <- res.open_ssd$SetInfo[res.open_ssd$SetInfo$SetID == gene.name, "SetIndex"]
  geno.mat <- Get_Genotypes_SSD(SSD_INFO = res.open_ssd,
                                Set_Index = index, is_ID = TRUE)
  #res.skato <- SKATBinary_Robust(geno.mat, obj = model, method = "SKATO")
  res.skato <- SKATBinary_Robust(geno.mat, obj = model, method = "SKATO", impute.method = impute.method)
  return(res.skato)
}
ALLcase.res.cohort1 <- ALLcase.SKATBinary_robust(pop= "ALLcasecontrol_cohort1", gene.name="EVC2")
ALLcase.res.cohort2 <- ALLcase.SKATBinary_robust(pop= "ALLcasecontrol_cohort2", gene.name="EVC2")







  

# names_case_cohort <- getDiseaseCohortNames(dataFile = pop)
# filepath="/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/SSD_results/"
# res.open_ssd <- Open_SSD(File.SSD = paste0(filepath, names_case_cohort$disease_name,"casecontrol_", names_case_cohort$cohort_name,"/", "File.SSD"),
#                          File.Info = paste0(filepath, names_case_cohort$disease_name,"casecontrol_", names_case_cohort$cohort_name,"/","File.Info"))
# 
# filepath.model <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/NullModels_OrderedPheno/"
# load(file=paste0(filepath.model, "nullModel_orderedPheno_", pop, ".RData"))
# phenotype_data <- res.phenotype_nullModel$phenotype_Data
# pheno_with_FIDIID <- res.phenotype_nullModel$pheno_with_FIDIID
# index <- res.open_ssd$SetInfo[res.open_ssd$SetInfo$SetID == gene.name, "SetIndex"]
# geno.mat <- Get_Genotypes_SSD(SSD_INFO = res.open_ssd,
#                               Set_Index = index, is_ID = TRUE)
# 
# 
# #--------------------------------T-cell SKATBinary_robust------------------------------------------
# control.IID <- pheno_with_FIDIID$IID[pheno_with_FIDIID$phenotype == 0]
# Tcell.casecontrol.IID <- c(case.Tcell.IID, control.IID)
# # phenotype data
# casecontrol.Tcell.dat <- phenotype_data[!is.na(match(pheno_with_FIDIID$IID, Tcell.casecontrol.IID)), ]
# # genotype data
# geno.mat.Tcell <- geno.mat[!is.na(match(pheno_with_FIDIID$IID, Tcell.casecontrol.IID)), ]
# 
# # refit-null model
# model <- SKAT_Null_Model(phenotype ~ .,
#                          data = casecontrol.Tcell.dat , out_type = "D")
# # flag <- apply(geno.mat.Tcell , 1, function(x) {any(is.na(x) | x == 9)})
# res.skato <- SKATBinary_Robust(geno.mat.Tcell, obj = model, method = "SKATO")
# #--------------------------------------------------------------------------------------------------
# 
# 
# #--------------------------------B-cell SKATBinary_robust------------------------------------------
# Bcell.casecontrol.IID <- c(case.Bcell.IID, control.IID)
# # phenotype data
# casecontrol.Bcell.dat <- phenotype_data[!is.na(match(pheno_with_FIDIID$IID, Bcell.casecontrol.IID)), ]
# # genotype data
# geno.mat.Bcell <- geno.mat[!is.na(match(pheno_with_FIDIID$IID, Bcell.casecontrol.IID)), ]
# 
# # refit-null model
# model <- SKAT_Null_Model(phenotype ~ .,
#                          data = casecontrol.Bcell.dat , out_type = "D")
# # flag <- apply(geno.mat.Tcell , 1, function(x) {any(is.na(x) | x == 9)})
# res.skato <- SKATBinary_Robust(geno.mat.Bcell, obj = model, method = "SKATO")
# #--------------------------------------------------------------------------------------------------
# 
# 
