
rm(list = ls())
library(MetaSKAT)
library(SKAT)

getDiseaseCohortNames <- function(dataFile="ALLcasecontrol_cohort1") {
  split.Full_name <- unlist(strsplit(dataFile, split="_"))
  cohort.name <- split.Full_name[2]
  disease.name <- "AMLMDS"
  if (paste(unlist(strsplit(split.Full_name[1], split =""))[1:3], collapse = "") == "ALL") {
    disease.name <- "ALL"
  }
  return(list("disease_name"=disease.name, "cohort_name"=cohort.name))
}

filepath.nullModel <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/subtypes_Tcell_Bcell/"
load(file=paste0(filepath.nullModel, paste0("ALL.EVC2.cohort1_nullModel_testRes_sampleList.RData")))
load(file=paste0(filepath.nullModel, paste0("ALL.EVC2.cohort2_nullModel_testRes_sampleList.RData")))




gen_MSSD_forSubtypes <- function(i, subtype="Tcell") {
  names_case_cohort <- getDiseaseCohortNames(dataFile = pop[i])
  cohort = names_case_cohort$cohort_name
  
  if (cohort == "cohort1") {
    subtype.res <- subtype.res.cohort1
  } else if (cohort == "cohort2") {
    subtype.res <- subtype.res.cohort2
  }
  
  if (subtype == "Tcell") {
    nullModel <- subtype.res$nullModel.Tcell
    n.sample <- dim(subtype.res$Tcell.casecontrol.sampleList)[1]
    sampleList <- subtype.res$Tcell.casecontrol.sampleList
  } else if (subtype == "Bcell") {
    nullModel <- subtype.res$nullModel.Bcell
    n.sample <- dim(subtype.res$Bcell.casecontrol.sampleList)[1]
    sampleList <- subtype.res$Bcell.casecontrol.sampleList
  }
  
  # check if mismatch
  filepath.bin=paste0("/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/subtypes_Tcell_Bcell/", subtype, "_Meta/")
  File.Fam = paste0(filepath.bin, "CLEAN.rare.", names_case_cohort$disease_name,"casecontrol_", names_case_cohort$cohort_name, ".fam")
  FAM <- Read_Plink_FAM(File.Fam)
  print(paste0("For pop ", pop[i], ", subtype= ", subtype, ", mismatch count is ", sum(sampleList$IID != FAM$IID)))
  
  filepath.bin=paste0("/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/subtypes_Tcell_Bcell/", subtype, "_Meta/")
  filepath.SetID="/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/SSD_results/"
  filepath.META_Files=paste0("/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/subtypes_Tcell_Bcell/", subtype, "_Meta/")
  res.gen_SSD_MetaFiles <- Generate_Meta_Files(obj = nullModel ,
                                               File.Bed = paste0(filepath.bin, "CLEAN.rare.", names_case_cohort$disease_name,"casecontrol_", 
                                                                 names_case_cohort$cohort_name, ".bed"),
                                               File.Bim = paste0(filepath.bin, "CLEAN.rare.", names_case_cohort$disease_name,"casecontrol_", 
                                                                 names_case_cohort$cohort_name, ".bim"),
                                               File.SetID = paste0(filepath.SetID, "File.SetID.autosome"), 
                                               
                                               File.MSSD = paste0(filepath.META_Files, names_case_cohort$disease_name,"casecontrol_", 
                                                                  names_case_cohort$cohort_name,  "/File.MSSD"),
                                               File.MInfo = paste0(filepath.META_Files, names_case_cohort$disease_name,"casecontrol_", 
                                                                   names_case_cohort$cohort_name, "/File.MInfo"),
                                               N.Sample = n.sample)
}


pop <- c("ALLcasecontrol_cohort1", "ALLcasecontrol_cohort2")
# cohort1, Tcell
gen_MSSD_forSubtypes(i=1, subtype = "Tcell")
# cohort1, Bcell
gen_MSSD_forSubtypes(i=1, subtype = "Bcell")
# cohort2, Tcell
gen_MSSD_forSubtypes(i=2, subtype = "Tcell")
# cohort2, Bcell
gen_MSSD_forSubtypes(i=2, subtype = "Bcell")




# T-cell cohort1 + cohort2
filepath.Tcell <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/subtypes_Tcell_Bcell/Tcell_Meta/"
filepath.Tcell.cohort1 <- paste0(filepath.Tcell, "ALLcasecontrol_cohort1")
filepath.Tcell.cohort2 <- paste0(filepath.Tcell, "ALLcasecontrol_cohort2")
File.MSSD.vec <- c(paste0(filepath.Tcell.cohort1, "/File.MSSD"),
                   paste0(filepath.Tcell.cohort2, "/File.MSSD"))
File.MInfo.vec <- c(paste0(filepath.Tcell.cohort1, "/File.MInfo"),
                    paste0(filepath.Tcell.cohort2, "/File.MInfo"))
res.open_MSSD <- Open_MSSD_File_2Read(File.MSSD.vec = File.MSSD.vec, 
                                      File.MInfo.vec = File.MInfo.vec)
Info_MAF_Missing_Score_Allele12.Tcell.cohort1  <- res.open_MSSD$EachInfo[[1]]$Info
Info_MAF_Missing_Score_Allele12.Tcell.cohort2  <- res.open_MSSD$EachInfo[[2]]$Info
EVC2.Tcell.cohort1 <- Info_MAF_Missing_Score_Allele12.Tcell.cohort1[Info_MAF_Missing_Score_Allele12.Tcell.cohort1$SetID=="EVC2",]
EVC2.Tcell.cohort2 <- Info_MAF_Missing_Score_Allele12.Tcell.cohort2[Info_MAF_Missing_Score_Allele12.Tcell.cohort2$SetID=="EVC2",]
dim(EVC2.Tcell.cohort1)
EVC2.Tcell.cohort1
EVC2.Tcell.cohort2


# B-cell cohort1 + cohort2
filepath.Bcell <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/subtypes_Tcell_Bcell/Bcell_Meta/"
filepath.Bcell.cohort1 <- paste0(filepath.Bcell, "ALLcasecontrol_cohort1")
filepath.Bcell.cohort2 <- paste0(filepath.Bcell, "ALLcasecontrol_cohort2")
File.MSSD.vec <- c(paste0(filepath.Bcell.cohort1, "/File.MSSD"),
                   paste0(filepath.Bcell.cohort2, "/File.MSSD"))
File.MInfo.vec <- c(paste0(filepath.Bcell.cohort1, "/File.MInfo"),
                    paste0(filepath.Bcell.cohort2, "/File.MInfo"))
res.open_MSSD <- Open_MSSD_File_2Read(File.MSSD.vec = File.MSSD.vec, 
                                      File.MInfo.vec = File.MInfo.vec)
Info_MAF_Missing_Score_Allele12.Bcell.cohort1  <- res.open_MSSD$EachInfo[[1]]$Info
Info_MAF_Missing_Score_Allele12.Bcell.cohort2  <- res.open_MSSD$EachInfo[[2]]$Info
EVC2.Bcell.cohort1 <- Info_MAF_Missing_Score_Allele12.Bcell.cohort1[Info_MAF_Missing_Score_Allele12.Bcell.cohort1$SetID=="EVC2",]
EVC2.Bcell.cohort2 <- Info_MAF_Missing_Score_Allele12.Bcell.cohort2[Info_MAF_Missing_Score_Allele12.Bcell.cohort2$SetID=="EVC2",]
EVC2.Bcell.cohort1
EVC2.Bcell.cohort2


subtype.res.cohort1$Bcell.res$param$n.marker.name

