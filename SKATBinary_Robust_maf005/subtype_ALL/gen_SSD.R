# This file is used to generate SSD files
rm(list=ls())
mypaths <- .libPaths()
mypaths <- c(mypaths, "/projects/rpci/qzhu/xinyu/Software/R/R-4.0.2")
.libPaths(mypaths)
# load library
library(skatMeta)
library(survivalSKAT)
library(SKAT)

#################################################################################
# # make File.SetID autosome 1-22 + sex chr X, Y
# temp.data <- snpInfo.functional[, colnames(snpInfo.functional)[c(5,1)]]
# head(temp.data)
# filepath="/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/SSD_results/"
# write.table(temp.data, file = paste0(filepath, "File.SetID"), 
#             row.names = FALSE, col.names = FALSE, 
#             quote = FALSE, sep=" ")


# # make File.SetID autosome 1-22
# snpInfo_autosome <- snpInfo.functional[(snpInfo.functional$Chr != 'X' & snpInfo.functional$Chr != 'Y'),]
# head(snpInfo_autosome)
# filepath="/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/SSD_results/"
# temp.data <- snpInfo_autosome[, colnames(snpInfo_autosome)[c(5,1)]]
# write.table(temp.data, file = paste0(filepath, "File.SetID.autosome"),
#             row.names = FALSE, col.names = FALSE,
#             quote = FALSE, sep=" ")
##################################################################################


# ##################################################################################
# # prepare for the manhattan plot
# # split data by "gene" 
# library(dplyr)
# head(snpInfo.functional)
# 
# gene_level.res <- snpInfo.functional %>% 
#                     group_by(Genes) %>% summarise(Chr=first(Chr), BP=median(MapInfo)) %>%
#                     arrange(Chr)
# 
# head(gene_level.res)
# filepath <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/SKATO_test_results/Chr_BP_Genes/"
# write.table(gene_level.res, file=paste0(filepath, "chr_BP_Genes.txt"),
#             row.names = FALSE, quote=FALSE)
# ##################################################################################





getDiseaseCohortNames <- function(dataFile="ALLcasecontrol_cohort1") {
  split.Full_name <- unlist(strsplit(dataFile, split="_"))
  cohort.name <- split.Full_name[2]
  disease.name <- "AMLMDS"
  if (paste(unlist(strsplit(split.Full_name[1], split =""))[1:3], collapse = "") == "ALL") {
    disease.name <- "ALL"
  }
  return(list("disease_name"=disease.name, "cohort_name"=cohort.name))
}

pop <- c("ALLcasecontrol_cohort1", "ALLcasecontrol_cohort2")

get_SSDplusInfo_fromCLEANdata <- function(i) {
  names_case_cohort <- getDiseaseCohortNames(dataFile = pop[i])

  
  filepath.SetID="/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/SSD_results/"
  # T-cell
  filepath.bin.Tcell="/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/cleaned_dataset/cleaned_rare.Tcell_maf005/"
  filepath.res.Tcell="/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/SSD_rare_maf005/Tcell/"
  # B-cell
  filepath.bin.Bcell="/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/cleaned_dataset/cleaned_rare.Bcell_maf005/"
  filepath.res.Bcell="/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/SSD_rare_maf005/Bcell/"
  
  # T-cell
  res.gen_SSD_setID <- Generate_SSD_SetID(File.Bed = paste0(filepath.bin.Tcell, "CLEAN.rare005.Tcell.", names_case_cohort$disease_name,"casecontrol_", 
                                                            names_case_cohort$cohort_name, ".bed"),
                                          File.Bim = paste0(filepath.bin.Tcell, "CLEAN.rare005.Tcell.", names_case_cohort$disease_name,"casecontrol_", 
                                                            names_case_cohort$cohort_name, ".bim"),
                                          File.Fam = paste0(filepath.bin.Tcell, "CLEAN.rare005.Tcell.", names_case_cohort$disease_name,"casecontrol_", 
                                                            names_case_cohort$cohort_name, ".fam"),
                                          File.SetID = paste0(filepath.SetID, "File.SetID.autosome"), 
                                          File.SSD = paste0(filepath.res.Tcell, names_case_cohort$disease_name,"casecontrol_", 
                                                            names_case_cohort$cohort_name,"/","File.SSD"),
                                          File.Info = paste0(filepath.res.Tcell, names_case_cohort$disease_name,"casecontrol_", 
                                                             names_case_cohort$cohort_name,"/","File.Info"))
  print("T-cell SSD files done!")
  
  # B-cell
  res.gen_SSD_setID <- Generate_SSD_SetID(File.Bed = paste0(filepath.bin.Bcell, "CLEAN.rare005.Bcell.", names_case_cohort$disease_name,"casecontrol_", 
                                                            names_case_cohort$cohort_name, ".bed"),
                                          File.Bim = paste0(filepath.bin.Bcell, "CLEAN.rare005.Bcell.", names_case_cohort$disease_name,"casecontrol_", 
                                                            names_case_cohort$cohort_name, ".bim"),
                                          File.Fam = paste0(filepath.bin.Bcell, "CLEAN.rare005.Bcell.", names_case_cohort$disease_name,"casecontrol_", 
                                                            names_case_cohort$cohort_name, ".fam"),
                                          File.SetID = paste0(filepath.SetID, "File.SetID.autosome"), 
                                          File.SSD = paste0(filepath.res.Bcell, names_case_cohort$disease_name,"casecontrol_", 
                                                            names_case_cohort$cohort_name,"/","File.SSD"),
                                          File.Info = paste0(filepath.res.Bcell, names_case_cohort$disease_name,"casecontrol_", 
                                                             names_case_cohort$cohort_name,"/","File.Info"))
  print("B-cell SSD files done!")
  
}




get_SSDplusInfo_fromCLEANdata(1)
get_SSDplusInfo_fromCLEANdata(2)



# till now, SSD file and Info file have been generated.
#####################################################################################################

# i=1
# names_case_cohort <- getDiseaseCohortNames(dataFile = pop[i])
# filepath.bin="/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/cleaned_dataset/cleaned_rare_variants/"
# File.Fam = paste0(filepath.bin, "CLEAN.rare.", names_case_cohort$disease_name,"casecontrol_", names_case_cohort$cohort_name, ".fam")
# FAM <- Read_Plink_FAM(File.Fam)
# 
# # load full phenotype data
# filepath.pheno <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/Penotype_Data/phenotypeData_without_Step7_2/"
# pheno.data <- read.table(paste0(filepath.pheno, "phenotype.data.", names_case_cohort$disease_name,
#                                 ".casecontrol.", names_case_cohort$cohort_name),
#                          header = T)
# # change order
# pheno.data <- pheno.data[match(FAM$IID, pheno.data$IID), ]





# get total number of sets in each pop, and divide these sets into several groups for batch computing
getDiseaseCohortNames <- function(dataFile="ALLcasecontrol_cohort1") {
  split.Full_name <- unlist(strsplit(dataFile, split="_"))
  cohort.name <- split.Full_name[2]
  disease.name <- "AMLMDS"
  if (paste(unlist(strsplit(split.Full_name[1], split =""))[1:3], collapse = "") == "ALL") {
    disease.name <- "ALL"
  }
  return(list("disease_name"=disease.name, "cohort_name"=cohort.name))
}

getNum_Sets <- function(i) {
  names_case_cohort <- getDiseaseCohortNames(dataFile = pop[i])
  filepath="/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/SSD_results/"
  res.open_ssd <- Open_SSD(File.SSD = paste0(filepath, names_case_cohort$disease_name,"casecontrol_", names_case_cohort$cohort_name,"/", "File.SSD"),
                           File.Info = paste0(filepath, names_case_cohort$disease_name,"casecontrol_", names_case_cohort$cohort_name,"/","File.Info"))
  
  filepath.out <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/SKATO_test_results/"
  write.table(res.open_ssd$nSets, file = paste0(filepath.out, pop[i], "/", "number_Sets.txt"),
              row.names = FALSE, col.names =FALSE, quote=FALSE)
  
  Close_SSD()
}
# getNum_Sets(1)
# getNum_Sets(2)
# getNum_Sets(3)
# getNum_Sets(4)


