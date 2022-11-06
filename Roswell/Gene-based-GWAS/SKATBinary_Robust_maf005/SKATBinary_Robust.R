#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# 1st arg: pop, e.g. ALLcasecontrol_cohort1
# 2nd & 3rd arg: start_index & end_index
pop = args[1]
start_index = args[2]
end_index = args[3]
# SKATBinary_Robust is for ALL cohort1 and 2


# this script receives 2 parameters: start_index and end_index of setID
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

names_case_cohort <- getDiseaseCohortNames(dataFile = pop)
filepath="/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/SSD_rare_maf005/"
res.open_ssd <- Open_SSD(File.SSD = paste0(filepath, names_case_cohort$disease_name,"casecontrol_", names_case_cohort$cohort_name,"/", "File.SSD"),
                         File.Info = paste0(filepath, names_case_cohort$disease_name,"casecontrol_", names_case_cohort$cohort_name,"/","File.Info"))
# res.open_ssd$nSets # total number of sets/Genes in this pop



############################################################################################
# pre-load phenotype data and null model
filepath.model <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/NullModels_OrderedPheno/"
load(file=paste0(filepath.model, "nullModel_orderedPheno_", pop, ".RData"))
phenotype_data <- res.phenotype_nullModel$phenotype_Data
null_model <- res.phenotype_nullModel$null_model
#pheno_with_FIDIID <- res.phenotype_nullModel$pheno_with_FIDIID
############################################################################################

# create lists to store results
start_time <- Sys.time()
test_Pvalue_list <- c()
set_list <- c()
testRes_RData_list <- list()

#start_index<-1501;end_index<-2000
iter <- 1
#start_index<-1;end_index<- 110

# loop from start_index to end_index
for (index in c(start_index:end_index)) {
  if (iter %% 20 == 0) {
    print(paste0("Job: ", iter))
  }
  
  # set/gene
  gene.name <- res.open_ssd$SetInfo[res.open_ssd$SetInfo$SetIndex == index, "SetID"]
  
  # load null_model and generate genotype matrix
  model <- null_model
  geno.mat <- Get_Genotypes_SSD(SSD_INFO = res.open_ssd,
                                Set_Index = index)
  
  ########################################################################
  # since MetaSkat does not provide the method which can remove missing
  # we have to use "fixed" impute method here
  
  # delete no-variants SNPs (columns)
  no_var_snps <- which(colSums(geno.mat) == 0)
  if (length(no_var_snps) > 0) {
    print("----------------------------------------------------------------------")
    print(paste0("Set index ", index, " gene_name ", gene.name, 
                 " has ", length(no_var_snps), " snp(s) that they have no variations after deleting missing samples!"))
    geno.mat <- geno.mat[, -no_var_snps]
    geno.mat <- as.matrix(geno.mat)
    print(paste0("Sex index ", index, " gene name ", gene.name,
                 " has removed columns with no variations"))
    print("----------------------------------------------------------------------")
  }

  ########################################################################
  # we will perform SKATO test only for genes with #n of snps >= 2
  if (dim(geno.mat)[2] < 2) {
    print(paste0("Sex index ", index, " gene name ", gene.name,
                 " has less than 2 variants! skip this gene."))
    iter <- iter+1
    next;
    
  }
  ########################################################################
  
  
  # SKATO test 
  # check orders of sample 
  res.skato <- SKATBinary_Robust(geno.mat, obj = model, method = "SKATO")
  
  # add this gene name into the list
  set_list <- c(set_list, gene.name)
  # add (detailed) test results into RData
  testRes_RData_list[[gene.name]] <- res.skato
  # add (gene.name and p-value) test results into txt
  test_Pvalue <- res.skato$p.value
  test_Pvalue_list <- c(test_Pvalue_list, test_Pvalue)
  
  iter <- iter + 1
}

Close_SSD()

print("For Loop is Done!")
end_time <- Sys.time()
print(paste0("cost time: ", end_time - start_time))

#---------------------------------------store results: Start!---------------------------------------#
# STORE results of SKAT_particular_Gene (RData) 
filepath.out <-"/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/SKATO_test_results_SKATBinary_Robust_BG_maf005/"
save(testRes_RData_list, file=paste0(filepath.out, pop, "/full_results/", 
                                     pop,"_full_testRes", start_index, "_", end_index,".RData"))

# store summary results and write into txt file
summary_output <- data.frame("SetID"=set_list,
                             "Pvalue"=test_Pvalue_list)
write.table(summary_output , paste0(filepath.out, pop, "/summary_results/", 
                                    "summary", start_index, "_", end_index, ".txt"),
            row.names = FALSE, quote = FALSE)
#---------------------------------------store results: Done!------------------------------------------#

print(paste0("Store results Job", start_index, "-", end_index, " is done!"))
