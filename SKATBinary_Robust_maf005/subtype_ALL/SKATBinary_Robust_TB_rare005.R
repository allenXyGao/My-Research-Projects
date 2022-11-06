#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# 1st arg: pop, e.g. ALLcasecontrol_cohort1
# 2nd & 3rd arg: start_index & end_index
pop = args[1]
start_index = args[2]
end_index = args[3]


library(SKAT)

# pop <- "ALLcasecontrol_cohort1"

getDiseaseCohortNames <- function(dataFile="ALLcasecontrol_cohort1") {
  split.Full_name <- unlist(strsplit(dataFile, split="_"))
  cohort.name <- split.Full_name[2]
  disease.name <- "AMLMDS"
  if (paste(unlist(strsplit(split.Full_name[1], split =""))[1:3], collapse = "") == "ALL") {
    disease.name <- "ALL"
  }
  return(list("disease_name"=disease.name, "cohort_name"=cohort.name))
}

names_case_cohort <- getDiseaseCohortNames(dataFile = pop)
filepath.Tcell <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/SSD_rare_maf005/Tcell/"
filepath.Bcell <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/SSD_rare_maf005/Bcell/"

res.open_ssd.Tcell <- Open_SSD(File.SSD = paste0(filepath.Tcell, names_case_cohort$disease_name,"casecontrol_", 
                                           names_case_cohort$cohort_name,"/", "File.SSD"),
                                File.Info = paste0(filepath.Tcell, names_case_cohort$disease_name,"casecontrol_", 
                                            names_case_cohort$cohort_name,"/","File.Info"))





############################################################################################
# pre-load phenotype data and null model
filepath.model <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/NullModels_OrderedPheno_Tcell_B_cell/"
load(file=paste0(filepath.model, "nullModel_orderedPheno_Tcell_Bcell.", pop,".RData"))
# load T cell
null_model.Tcell <- res.phenotype_nullModel$null_model_T.cell
# load B cell
null_model.Bcell <- res.phenotype_nullModel$null_model_B.cell
############################################################################################




# create lists to store results
test_Pvalue_list.Tcell <- c()
set_list.Tcell <- c()

iter <- 1
# start_index<-1;end_index<-100

# loop from start_index to end_index
for (index in c(start_index:end_index)) {
  if (iter %% 20 == 0) {
    print(paste0("Job: ", iter))
  }
  
  # set/gene
  gene.name <- res.open_ssd.Tcell$SetInfo[res.open_ssd.Tcell$SetInfo$SetIndex == index, "SetID"]
  geno.mat.Tcell <- Get_Genotypes_SSD(SSD_INFO = res.open_ssd.Tcell,
                                      Set_Index = index, is_ID = TRUE)
  
  
  # delete no-variants SNPs (columns)
  no_var_snps <- which(colSums(geno.mat.Tcell) == 0)
  if (length(no_var_snps) > 0) {
    print("------------------------T-cell-------------------------------------")
    print(paste0("Set index ", index, " gene_name ", gene.name, 
                 " has ", length(no_var_snps), " snp(s) that they have no variations after deleting missing samples!"))
    geno.mat.Tcell <- geno.mat.Tcell[, -no_var_snps]
    geno.mat.Tcell <- as.matrix(geno.mat.Tcell)
    print(paste0("Sex index ", index, " gene name ", gene.name,
                 " has removed columns with no variations"))
    print("------------------------T-cell-----------------------------------------")
  }

  
  ########################################################################
  # we will perform SKATO test only for genes with #n of snps >= 2
  if (dim(geno.mat.Tcell)[2] < 2) {
    print(paste0("T-cell, Sex index ", index, " gene name ", gene.name,
                 " has less than 2 variants! skip this gene."))
    p.val.Tcell <- NA
  } else {
    res.skato.Tcell <- SKATBinary_Robust(geno.mat.Tcell, obj = null_model.Tcell, method = "SKATO")
    p.val.Tcell <- res.skato.Tcell$p.value
  }
  test_Pvalue_list.Tcell <- c(test_Pvalue_list.Tcell, p.val.Tcell)
  ########################################################################
  
  # add this gene name into the list
  set_list.Tcell <- c(set_list.Tcell, gene.name)
  iter <- iter + 1

}

Close_SSD()
print("For T-cell Loop is Done!")



#-----------------------------------------------------------------------------------------------------#

res.open_ssd.Bcell <- Open_SSD(File.SSD = paste0(filepath.Bcell, names_case_cohort$disease_name,"casecontrol_", 
                                                 names_case_cohort$cohort_name,"/", "File.SSD"),
                               File.Info = paste0(filepath.Bcell, names_case_cohort$disease_name,"casecontrol_", 
                                                  names_case_cohort$cohort_name,"/","File.Info"))
test_Pvalue_list.Bcell <- c()
set_list.Bcell <- c()

iter <- 1
# start_index<-1;end_index<-100

# loop from start_index to end_index
for (index in c(start_index:end_index)) {
  if (iter %% 20 == 0) {
    print(paste0("Job: ", iter))
  }
  
  # set/gene
  gene.name <- res.open_ssd.Tcell$SetInfo[res.open_ssd.Tcell$SetInfo$SetIndex == index, "SetID"]
  geno.mat.Bcell <- Get_Genotypes_SSD(SSD_INFO = res.open_ssd.Bcell,
                                      Set_Index = index, is_ID = TRUE)
  
  no_var_snps <- which(colSums(geno.mat.Bcell) == 0)
  if (length(no_var_snps) > 0) {
    print("------------------------B-cell-------------------------------------")
    print(paste0("Set index ", index, " gene_name ", gene.name, 
                 " has ", length(no_var_snps), " snp(s) that they have no variations after deleting missing samples!"))
    geno.mat.Bcell <- geno.mat.Bcell[, -no_var_snps]
    geno.mat.Bcell <- as.matrix(geno.mat.Bcell)
    print(paste0("Sex index ", index, " gene name ", gene.name,
                 " has removed columns with no variations"))
    print("------------------------B-cell-----------------------------------------")
  }
  
  ########################################################################
  # we will perform SKATO test only for genes with #n of snps >= 2
  if (dim(geno.mat.Bcell)[2] < 2) {
    print(paste0("B-cell, Sex index ", index, " gene name ", gene.name,
                 " has less than 2 variants! skip this gene."))
    p.val.Bcell <- NA
  } else {
    res.skato.Bcell <- SKATBinary_Robust(geno.mat.Bcell, obj = null_model.Bcell, method = "SKATO")
    p.val.Bcell <- res.skato.Bcell$p.value
  }
  test_Pvalue_list.Bcell <- c(test_Pvalue_list.Bcell, p.val.Bcell)
  ########################################################################
  

  # add this gene name into the list
  set_list.Bcell <- c(set_list.Bcell, gene.name)
  iter <- iter + 1
  
}

Close_SSD()
print("For B-cell Loop is Done!")


#---------------------------------------store results: Start!---------------------------------------#
# store summary results and write into txt file

filepath.out <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/SKATO_test_results_SKATBinary_Robust_BG_maf005/Tcell/"
summary_output.Tcell <- data.frame("SetID"=set_list.Tcell,
                             "Pvalue"=test_Pvalue_list.Tcell)
write.table(summary_output.Tcell, paste0(filepath.out, pop, "/summary_results/", 
                                   "summary", start_index, "_", end_index, ".txt"),
            row.names = FALSE, quote = FALSE)


filepath.out <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/SKATO_test_results_SKATBinary_Robust_BG_maf005/Bcell/"
summary_output.Bcell <- data.frame("SetID"=set_list.Bcell,
                             "Pvalue"=test_Pvalue_list.Bcell)
write.table(summary_output.Bcell, paste0(filepath.out, pop, "/summary_results/", 
                                   "summary", start_index, "_", end_index, ".txt"),
            row.names = FALSE, quote = FALSE)
  
  
#---------------------------------------store results: Done!------------------------------------------#


print(paste0("Store results Job", start_index, "-", end_index, " is done!"))









