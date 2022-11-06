
# please the grab "amstatpr", "mdstatpr" columns from the meta-pheno file 
# /projects/rpci/qzhu/lsuchest/lsuchest/DBMT_PhenoData/DBMT_PhenoData_long_allVar_20180424.txt 
# for the AML patients carried DNMT3A variants 
# (the red ones in your variant plots for AMLMDS analysis).
library(SKAT)

filepath <- "/projects/rpci/qzhu/lsuchest/lsuchest/DBMT_PhenoData/"
d <- read.table(file = paste0(filepath, "DBMT_PhenoData_long_allVar_20180424.txt"), sep="\t", header = T) 
names.extract <- c("Exome.IID","Exome.FID", "amstatpr", "mdstatpr", "cohort", "disease")
d.extract <- d[, names.extract]
head(d.extract)
out <- strsplit(d.extract$Exome.IID, "([A-Z]-)", perl = TRUE)
d.extract$ID <- gsub("[-_]", "", do.call(rbind, out)[,2])


# find ID of AML patients in cohort1 carrying variants in DNMT3A  
getDiseaseCohortNames <- function(dataFile="ALLcasecontrol_cohort1") {
  split.Full_name <- unlist(strsplit(dataFile, split="_"))
  cohort.name <- split.Full_name[2]
  disease.name <- "AMLMDS"
  if (paste(unlist(strsplit(split.Full_name[1], split =""))[1:3], collapse = "") == "ALL") {
    disease.name <- "ALL"
  }
  return(list("disease_name"=disease.name, "cohort_name"=cohort.name))
}
#######################################################################################################################
pop <- "AMLMDScasecontrol_cohort1"
gene.name <- "DNMT3A"
#######################################################################################################################
names_case_cohort <- getDiseaseCohortNames(dataFile = pop)
filepath="/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/SSD_results/"
res.open_ssd <- Open_SSD(File.SSD = paste0(filepath, names_case_cohort$disease_name,"casecontrol_", names_case_cohort$cohort_name,"/", "File.SSD"),
                         File.Info = paste0(filepath, names_case_cohort$disease_name,"casecontrol_", names_case_cohort$cohort_name,"/","File.Info"))

filepath.model <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/NullModels_OrderedPheno/"
load(file=paste0(filepath.model, "nullModel_orderedPheno_", pop, ".RData"))
phenotype_data <- res.phenotype_nullModel$phenotype_Data
pheno_with_FIDIID <- res.phenotype_nullModel$pheno_with_FIDIID

case.control_status <- data.frame("ID"=pheno_with_FIDIID$IID, "status"=pheno_with_FIDIID$phenotype)
out <- strsplit(case.control_status$ID, "([A-Z]-)", perl = TRUE)
case.control_status$ID <- gsub("[-_]", "", do.call(rbind, out)[,2])

index <- res.open_ssd$SetInfo[res.open_ssd$SetInfo$SetID == gene.name, "SetIndex"]
geno.mat <- Get_Genotypes_SSD(SSD_INFO = res.open_ssd,
                              Set_Index = index, is_ID = TRUE)
carriers <- apply(geno.mat, 2, function(x) which(x!=0 & x!=9, arr.ind = T))
inds <- unique(unlist(carriers))
inds.ordered <- inds[order(inds)]
# case/control
ind.case <- match(inds.ordered, which(case.control_status$status == 1))
ind.case = inds.ordered[!is.na(ind.case)]


patients.IID.AMLMDS.carried.DNMT3A <- case.control_status$ID[ind.case]
res <- d.extract[!is.na(match(d.extract$ID, patients.IID.AMLMDS.carried.DNMT3A)), ]
res$ID <- as.character(res$ID)
res
# store results
filepath <- "/user/xinyugao/SKATO_Project/grab_amstatpr_mdstatpr_columns/"
write.table(res, file = paste0(filepath, pop, "_", gene.name, "_info.amstatpr_mdstatpr.txt")
            ,row.names = FALSE, quote = TRUE)

# store results Done!
#####################################################################################################

# read stored results
# cohort1
filepath <- "/user/xinyugao/SKATO_Project/grab_amstatpr_mdstatpr_columns/"
pop <- "AMLMDScasecontrol_cohort1"
gene.name <- "DNMT3A"
info.amstatpr_mdstatpr.DNMT3A.AMLMDS.cohort1 <- read.table(file = paste0(filepath, pop, "_", gene.name, "_info.amstatpr_mdstatpr.txt"), 
                                                           header = TRUE, 
                                                           colClasses = c(rep("character",4), "numeric", rep("character",2)))
# specify colClasses in case read.table truncate the leading zeros!!
info.amstatpr_mdstatpr.DNMT3A.AMLMDS.cohort1
dim(info.amstatpr_mdstatpr.DNMT3A.AMLMDS.cohort1)
# cohort2
pop <- "AMLMDScasecontrol_cohort2"
gene.name <- "DNMT3A"
info.amstatpr_mdstatpr.DNMT3A.AMLMDS.cohort2 <- read.table(file = paste0(filepath, pop, "_", gene.name, "_info.amstatpr_mdstatpr.txt"), 
                                                           header = TRUE,
                                                           colClasses = c(rep("character",4), "numeric", rep("character",2)))
info.amstatpr_mdstatpr.DNMT3A.AMLMDS.cohort2
dim(info.amstatpr_mdstatpr.DNMT3A.AMLMDS.cohort2)


#####################################################################################################

# simplify status -> CR Relapse NA
# cohort1
info.cohort1.merged <- rep(NA, dim(info.amstatpr_mdstatpr.DNMT3A.AMLMDS.cohort1)[1])
for (i in 1:dim(info.amstatpr_mdstatpr.DNMT3A.AMLMDS.cohort1)[1]) {
  if (!is.na(info.amstatpr_mdstatpr.DNMT3A.AMLMDS.cohort1$amstatpr[i])) {
    info.cohort1.merged[i] <- info.amstatpr_mdstatpr.DNMT3A.AMLMDS.cohort1$amstatpr[i]
  }
  if (!is.na(info.amstatpr_mdstatpr.DNMT3A.AMLMDS.cohort1$mdstatpr[i])) {
    info.cohort1.merged[i] <- info.amstatpr_mdstatpr.DNMT3A.AMLMDS.cohort1$mdstatpr[i]
  }
}
table(info.cohort1.merged)
info.amstatpr_mdstatpr.DNMT3A.AMLMDS.cohort1$status <- info.cohort1.merged
table(info.amstatpr_mdstatpr.DNMT3A.AMLMDS.cohort1$status)

info.amstatpr_mdstatpr.DNMT3A.AMLMDS.cohort1$status[info.amstatpr_mdstatpr.DNMT3A.AMLMDS.cohort1$status == "1st Complete remission" |
                                                      info.amstatpr_mdstatpr.DNMT3A.AMLMDS.cohort1$status == "2nd CR" |
                                                      info.amstatpr_mdstatpr.DNMT3A.AMLMDS.cohort1$status == "CR"] <- "CR"
info.amstatpr_mdstatpr.DNMT3A.AMLMDS.cohort1$status[info.amstatpr_mdstatpr.DNMT3A.AMLMDS.cohort1$status == "1st Relapse" |
                                                      info.amstatpr_mdstatpr.DNMT3A.AMLMDS.cohort1$status == "2nd relapse" |
                                                      info.amstatpr_mdstatpr.DNMT3A.AMLMDS.cohort1$status == "No treatment" |
                                                      info.amstatpr_mdstatpr.DNMT3A.AMLMDS.cohort1$status =="Primary induction failure"] <- "Relapse"
info.status.cohort1 <- info.amstatpr_mdstatpr.DNMT3A.AMLMDS.cohort1[, c("ID", "status")]
# store simplified status
filepath <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/variant_plots/"
write.table(info.status.cohort1 , file = paste0(filepath, gene.name, "_info.status_cohort1.txt")
            ,row.names = FALSE)


# cohort2
info.cohort2.merged <- rep(NA, dim(info.amstatpr_mdstatpr.DNMT3A.AMLMDS.cohort2)[1])
for (i in 1:dim(info.amstatpr_mdstatpr.DNMT3A.AMLMDS.cohort2)[1]) {
  if (!is.na(info.amstatpr_mdstatpr.DNMT3A.AMLMDS.cohort2$amstatpr[i])) {
    info.cohort2.merged[i] <- info.amstatpr_mdstatpr.DNMT3A.AMLMDS.cohort2$amstatpr[i]
  }
  if (!is.na(info.amstatpr_mdstatpr.DNMT3A.AMLMDS.cohort2$mdstatpr[i])) {
    info.cohort2.merged[i] <- info.amstatpr_mdstatpr.DNMT3A.AMLMDS.cohort2$mdstatpr[i]
  }
}
table(info.cohort2.merged)
info.amstatpr_mdstatpr.DNMT3A.AMLMDS.cohort2$status <- info.cohort2.merged
table(info.amstatpr_mdstatpr.DNMT3A.AMLMDS.cohort2$status)

info.amstatpr_mdstatpr.DNMT3A.AMLMDS.cohort2$status[info.amstatpr_mdstatpr.DNMT3A.AMLMDS.cohort2$status == "1st Complete remission" |
                                                      info.amstatpr_mdstatpr.DNMT3A.AMLMDS.cohort2$status == "2nd CR" |
                                                      info.amstatpr_mdstatpr.DNMT3A.AMLMDS.cohort2$status == "CR"] <- "CR"
info.amstatpr_mdstatpr.DNMT3A.AMLMDS.cohort2$status[info.amstatpr_mdstatpr.DNMT3A.AMLMDS.cohort2$status == "1st Relapse" |
                                                      info.amstatpr_mdstatpr.DNMT3A.AMLMDS.cohort2$status == "2nd relapse" |
                                                      info.amstatpr_mdstatpr.DNMT3A.AMLMDS.cohort2$status == "No treatment" |
                                                      info.amstatpr_mdstatpr.DNMT3A.AMLMDS.cohort2$status =="Primary induction failure"] <- "Relapse"
info.status.cohort2 <- info.amstatpr_mdstatpr.DNMT3A.AMLMDS.cohort2[, c("ID", "status")]
table(info.status.cohort2$status)
filepath <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/variant_plots/"
write.table(info.status.cohort2 , file = paste0(filepath, gene.name, "_info.status_cohort2.txt")
            ,row.names = FALSE)
# store simplified status done!









