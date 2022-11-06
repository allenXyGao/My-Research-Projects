
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
pop <- "AMLMDScasecontrol_cohort2"
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
ind.control <- match(inds.ordered, which(case.control_status$status == 0))
ind.control = inds.ordered[!is.na(ind.control)]

# control group: carriers
control.IID.AMLMDS.carried.DNMT3A <- case.control_status$ID[ind.control]

patients.IID.AMLMDS.carried.DNMT3A <- case.control_status$ID[ind.case]
res <- d.extract[!is.na(match(d.extract$ID, patients.IID.AMLMDS.carried.DNMT3A)), ]
res$ID <- as.character(res$ID)
res

# patients not carriers of SNPs within DNMT3A
case.ID.not.carriers.DNMT3A <- case.control_status$ID[! case.control_status$ID %in% res$ID & case.control_status$status == 1 ]
length(case.ID.not.carriers.DNMT3A )
res.case.not.carriers <- d.extract[!is.na(match(d.extract$ID, case.ID.not.carriers.DNMT3A)), ]

info.merged <- rep(NA, dim(res.case.not.carriers)[1])
for (i in 1:dim(res.case.not.carriers)[1]) {
  if (!is.na(res.case.not.carriers$amstatpr[i])) {
    info.merged[i] <- res.case.not.carriers$amstatpr[i]
  }
  if (!is.na(res.case.not.carriers$mdstatpr[i])) {
    info.merged[i] <- res.case.not.carriers$mdstatpr[i]
  }
}
table(info.merged)
sum(is.na(info.merged))


## cohort1
# 1st Complete remission        1st Relapse                     2nd CR                       2nd relapse 
# 528                             131                             254                              51 
# 3rd CR                          CR                          ge 3rd CR                       ge 4th CR 
# 20                              11                               6                               1 
# Hematologic improvement (HI)      No response/stable disease         No treatment       Primary induction failure 
# 1                                       18                               8                             148 
# Relapse from complete remission   NA 
# 1                                 303

## cohort2
# 1st Complete remission      1st Relapse                       2nd CR                     2nd relapse 
# 180                              22                              57                               6 
# CR                       ge 3rd CR             Hematologic improvement (HI)      No response/stable disease 
# 18                               5                              28                              61 
# No treatment       Primary induction failure                Prog from HI        Relapse from complete remission 
# 26                              39                               8                               1 


res.case.not.carriers$status <- info.merged
res.case.not.carriers$status[res.case.not.carriers$status == "1st Complete remission" |
                               res.case.not.carriers$status == "2nd CR" | res.case.not.carriers$status == "3rd CR" |
                               res.case.not.carriers$status == "CR" | res.case.not.carriers$status == "ge 3rd CR" |
                               res.case.not.carriers$status == "ge 4th CR" ] <- "CR"

res.case.not.carriers$status[res.case.not.carriers$status == "1st Relapse" |
                               res.case.not.carriers$status == "2nd relapse" |
                               res.case.not.carriers$status == "No treatment" |
                               res.case.not.carriers$status =="Primary induction failure" |
                               res.case.not.carriers$status == "Relapse from complete remission" |
                               res.case.not.carriers$status == "No response/stable disease" |
                               res.case.not.carriers$status == "Hematologic improvement (HI)" |
                               res.case.not.carriers$status == "Prog from HI"] <- "Relapse"

table(res.case.not.carriers$status)

sum(is.na(res.case.not.carriers$status))



