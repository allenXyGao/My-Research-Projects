rm(list=ls())
library(dplyr)
library(tidyr)
library(SKAT)

# generate another version of the variant plots for ALL case-control where you color the cases 
# based on the "Lineage" column of file /projects/rpci/qzhu/lsuchest/lsuchest/DBMT_PhenoData/ALL_subtype_pheno_long.txt?

filepath <- "/projects/rpci/qzhu/lsuchest/lsuchest/DBMT_PhenoData/"
d <- read.table(file = paste0(filepath, "ALL_subtype_pheno_long.txt"), sep="\t", header = T) 
dim(d)
head(d)

names.select <- c("IID", "updated.EA_IID", "Lineage")
d.select <- d[, names.select]
head(d.select, 20)
table(d.select$Lineage)

# cohort1,2 subtype might be different
# Mature B-cell   precursor B-cell      T-cell (associated with ALL?)
#      14              822              146   NA
#      orange        xx                xx    0
#  order 1. sample by subtype when drawing
#        2. variants carried by samples

# out <- strsplit(d.select$IID, "([A-Z]-)", perl = TRUE)
# d.select$ID <- gsub("[-_]", "", do.call(rbind, out)[,2])
# d.select

#####################################################################################################################
# cohort.name <- "cohort1"
# disease.name <- "ALL"
# /projects/rpci/qzhu/lsuchest/lsuchest/DBMT_PhenoData/DBMT_PhenoData_long_allVar_20180424.txt
filepath <- "/projects/rpci/qzhu/lsuchest/lsuchest/DBMT_PhenoData/"
meta_phenotype <- read.table(file=paste(filepath, "DBMT_PhenoData_long_allVar_20180424.txt", sep=""), sep="\t", header = T)
colnames(meta_phenotype)
# ALL cohort 1
meta_phenotype.ALL.cohort1 <- meta_phenotype[meta_phenotype$Exome.recip.pop == "EA_cohort1" & meta_phenotype$disease == "ALL", 
                                             c("Exome.IID", "IID", "updated.EA_IID")]
head(meta_phenotype.ALL.cohort1, 20)
#ALL cohort 2
meta_phenotype.ALL.cohort2 <- meta_phenotype[meta_phenotype$Exome.recip.pop == "EA_cohort2" & meta_phenotype$disease == "ALL", 
                                             c("Exome.IID", "IID", "updated.EA_IID")]
head(meta_phenotype.ALL.cohort2, 20)

# meta_phenotype.ALL.cohort1
#################################################################################################################################


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
pop <- "ALLcasecontrol_cohort1"
gene.name <- "EVC2"
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
# out <- strsplit(case.control_status$ID, "([A-Z]-)", perl = TRUE)
# case.control_status$ID <- gsub("[-_]", "", do.call(rbind, out)[,2])
index <- res.open_ssd$SetInfo[res.open_ssd$SetInfo$SetID == gene.name, "SetIndex"]
geno.mat <- Get_Genotypes_SSD(SSD_INFO = res.open_ssd,
                              Set_Index = index, is_ID = TRUE)
carriers <- apply(geno.mat, 2, function(x) which(x!=0 & x!=9, arr.ind = T))
inds <- unique(unlist(carriers))
inds.ordered <- inds[order(inds)]
# case/control
ind.case <- match(inds.ordered, which(case.control_status$status == 1))
ind.case = inds.ordered[!is.na(ind.case)]
# ID of patients in cohort1 carrying EVC2
patients.IID.ALL.cohort1.carried.EVC2 <- case.control_status$ID[ind.case] 
# match
meta_phenotype.ALL.EVC2.cohort1 <- meta_phenotype.ALL.cohort1[!is.na(match(meta_phenotype.ALL.cohort1$Exome.IID, patients.IID.ALL.cohort1.carried.EVC2)),]
meta_phenotype.ALL.EVC2.cohort1
#-----------------------------------------------------------------------------------------------------------------------------------------#
# cohort2
pop2 <- "ALLcasecontrol_cohort2"
gene.name <- "EVC2"
names_case_cohort <- getDiseaseCohortNames(dataFile = pop2)
filepath="/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/SSD_results/"
res.open_ssd <- Open_SSD(File.SSD = paste0(filepath, names_case_cohort$disease_name,"casecontrol_", names_case_cohort$cohort_name,"/", "File.SSD"),
                         File.Info = paste0(filepath, names_case_cohort$disease_name,"casecontrol_", names_case_cohort$cohort_name,"/","File.Info"))

filepath.model <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/NullModels_OrderedPheno/"
load(file=paste0(filepath.model, "nullModel_orderedPheno_", pop2, ".RData"))
phenotype_data <- res.phenotype_nullModel$phenotype_Data
pheno_with_FIDIID <- res.phenotype_nullModel$pheno_with_FIDIID

case.control_status <- data.frame("ID"=pheno_with_FIDIID$IID, "status"=pheno_with_FIDIID$phenotype)
# out <- strsplit(case.control_status$ID, "([A-Z]-)", perl = TRUE)
# case.control_status$ID <- gsub("[-_]", "", do.call(rbind, out)[,2])
index <- res.open_ssd$SetInfo[res.open_ssd$SetInfo$SetID == gene.name, "SetIndex"]
geno.mat <- Get_Genotypes_SSD(SSD_INFO = res.open_ssd,
                              Set_Index = index, is_ID = TRUE)
carriers <- apply(geno.mat, 2, function(x) which(x!=0 & x!=9, arr.ind = T))
inds <- unique(unlist(carriers))
inds.ordered <- inds[order(inds)]
# case/control
ind.case <- match(inds.ordered, which(case.control_status$status == 1))
ind.case = inds.ordered[!is.na(ind.case)]
# ID of patients in cohort2 carrying EVC2
patients.IID.ALL.cohort2.carried.EVC2 <- case.control_status$ID[ind.case] 
# match
meta_phenotype.ALL.EVC2.cohort2 <- meta_phenotype.ALL.cohort2[!is.na(match(meta_phenotype.ALL.cohort2$Exome.IID, patients.IID.ALL.cohort2.carried.EVC2)),]
meta_phenotype.ALL.EVC2.cohort2

###############################################################################################################################################

# results
# cohort1
ALL.cohort1.EVC2.subtype <- meta_phenotype.ALL.EVC2.cohort1 %>% left_join(d.select, by=c("IID"))
ALL.cohort1.EVC2.subtype.res <- ALL.cohort1.EVC2.subtype[, c("Exome.IID", "Lineage")]

sum(is.na(ALL.cohort1.EVC2.subtype.res$Lineage))
table(ALL.cohort1.EVC2.subtype.res$Lineage)
# cohort2
ALL.cohort2.EVC2.subtype <- meta_phenotype.ALL.EVC2.cohort2 %>% left_join(d.select, by=c("IID"))
ALL.cohort2.EVC2.subtype.res <- ALL.cohort2.EVC2.subtype[, c("Exome.IID", "Lineage")]

# store results
file.out <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/variant_plots/"
write.table(ALL.cohort1.EVC2.subtype.res, file = paste0(file.out, "EVC2_Cohort1_Info.Lineage.txt")
            ,row.names = FALSE, quote = TRUE)
write.table(ALL.cohort2.EVC2.subtype.res, file = paste0(file.out, "EVC2_Cohort2_Info.Lineage.txt")
            ,row.names = FALSE, quote = TRUE)


###################################################################################################################################################
# store results done!


###################################################################################################################################################
# draw variant plots

file.out <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/variant_plots/"
info.status.cohort1 <- read.table(file = paste0(file.out, "EVC2_Cohort1_Info.Lineage.txt"), 
                                  header = TRUE)
info.status.cohort2 <- read.table(file = paste0(file.out, "EVC2_Cohort2_Info.Lineage.txt"), 
                                  header = TRUE)

pop <- "ALLcasecontrol_cohort1"
gene.name <- "EVC2"

names_case_cohort <- getDiseaseCohortNames(dataFile = pop)
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
snps <- colnames(geno.mat)
case.control_status <- data.frame("ID"=pheno_with_FIDIID$IID, "status"=pheno_with_FIDIID$phenotype)
# out <- strsplit(case.control_status$ID, "([A-Z]-)", perl = TRUE)
# case.control_status$ID <- gsub("[-_]", "", do.call(rbind, out)[,2])
head(case.control_status) # 1: case; 0: control

carriers <- apply(geno.mat, 2, function(x) which(x!=0 & x!=9, arr.ind = T))
inds <- unique(unlist(carriers))
inds.ordered <- inds[order(inds)]
# case/control
ind.case <- match(inds.ordered, which(case.control_status$status == 1))
ind.case = inds.ordered[!is.na(ind.case)]
ind.control <- match(inds.ordered, which(case.control_status$status == 0))
ind.control = inds.ordered[!is.na(ind.control)]
case.control.col <- c(rep("red", length(ind.case)), rep("blue", length(ind.control)))

n_snps <- length(snps); ninds <- length(inds)
col.mat <- matrix("grey", nrow = n_snps, ncol = ninds)

change_col <- function(col.mat, info.status) {
  for (i in 1:n_snps) {
    print(paste0("current snp is: ", snps[i]))
    #
    if (is.na(carriers[snps[i]])) {
      cur_carriers <- carriers
    } else {
      cur_carriers <- carriers[[i]]
    }
    
    for (j in 1:length(cur_carriers)) {
      matched_sample <- which(!is.na(match(inds.ordered, cur_carriers[j])))
      if (matched_sample <= length(ind.case)) {
        # further check if status is remission or relapse or NA
        ID <- case.control_status$ID[inds.ordered[matched_sample]]
        Lineage <- info.status$Lineage[info.status$Exome.IID == ID]
        if (is.na(Lineage)) {
          col.mat[i, matched_sample] = "red"
        }
        else if (Lineage == "Mature B-cell") {
          col.mat[i, matched_sample] = "plum"
        }
        else if (Lineage == "precursor B-cell") {
          col.mat[i, matched_sample] = "orange"
        }
        else if (Lineage == "T-cell") {
          col.mat[i, matched_sample] = "hotpink"
        }
      }
      else {
        col.mat[i, matched_sample] = "blue"
      }
    }
    print(paste0("current snp is: ", snps[i], " done!"))
  }
  return(col.mat)
}

col.mat <- change_col(col.mat = col.mat, info.status = info.status.cohort1)
#col.mat <- change_col(col.mat = col.mat, info.status = info.status.cohort2)
col.mat <- as.matrix(col.mat, nrow = n_snps, ncol = ninds)

# next reorder rows and columns
if (dim(col.mat)[1] == 1) {
  col.mat <- col.mat
} else {
  row.sums.raw <- rowSums(col.mat != "grey")
  row.order.info <- sort(row.sums.raw, decreasing = TRUE,index.return=TRUE)
  col.mat <- col.mat[row.order.info$ix, ] # update col.mat
  snps <- snps[row.order.info$ix]# update snp names
}


#### reorder columns
# each column can have only one color -> identify the category of the first non-grey element  
carriers_ID <- case.control_status$ID[inds.ordered]
column.status.info <- rep("grey", ninds)
for (j in 1:ninds) {
  match.info <- col.mat[,j][col.mat[, j] != "grey"]
  print(paste0("j= ", j, "match.info= ", match.info))
  if (is.na(match.info)) {next;}
  column.status.info[j] <- match.info[1]
}
#column.status.info
column.status.info <- factor(column.status.info, levels = c("orange", "hotpink", "red", "blue"))
column.order.index <- order(column.status.info)
column.status.info[column.order.index]

col.mat <- col.mat[, column.order.index] # update col.mat
carriers_ID <- carriers_ID[column.order.index] # update carriers ID
out <- strsplit(carriers_ID, "([A-Z]-)", perl = TRUE)
carriers_ID <- gsub("[-_]", "", do.call(rbind, out)[,2])
col.mat <- matrix(col.mat, nrow = n_snps, ncol = ninds)


setwd("/user/xinyugao/SKATO_Project/plots_carriers_variants/")
pdf(paste0("subtype_summary_", pop, "_",gene.name, ".pdf"))
xleft = 0:(ninds-1)
xright = xleft+1
vertical.len <- 1
plot(0, 0, type = "n", xlim = c(0, ninds), ylim = c(0, n_snps+1), 
     axes = FALSE, xlab = "", ylab = "" )
rect(xleft, n_snps/(1/vertical.len), xright, n_snps/(1/vertical.len)+vertical.len, border = "light gray" , col=case.control.col)
axis(2, at=n_snps/(1/vertical.len)+vertical.len/2, labels="case/control status", tick=F, lty=6, 
     pos=0.3, las=2, cex.axis=0.6)
axis(3, at=xleft+0.5, labels=carriers_ID, tick=F , pos=n_snps/(1/vertical.len)+vertical.len, 
     las=2, cex.axis=0.4, xpd = TRUE)
start = n_snps / (1/vertical.len)
ybottom <- start - seq(vertical.len, n_snps / (1/vertical.len), by = vertical.len)
ytop <- ybottom + vertical.len

for (j in 1:n_snps) {
  rect(xleft, ybottom[j], xright, ytop[j], border="light grey", col=col.mat[j, ])
  axis(2, at=ybottom[j]+vertical.len/2, labels = snps[j], tick=F, lty=6, pos=0, las=2, cex.axis=0.6)
  axis(4, at=ybottom[j]+vertical.len/2, labels=sum(col.mat[j, ]!="grey"), 
       tick=F, lty=6, pos=ninds-0.2, las=2, cex.axis=0.7)
}
axis(1, at=xleft+0.5, labels=colSums(col.mat != "grey"), tick=F, pos=0.1, las=2, cex.axis=0.5, xpd = TRUE)

# legend(x=0, y=-1.5, legend=c("precursor B-cell", "T-cell", "case (unknown)"), 
#        pt.bg=c("orange", "hotpink", "red"), pch=22, cex=0.6, xpd = TRUE)
legend(x=0, y=-0.8, legend=c("precursor B-cell"), 
       pt.bg=c("orange"), pch=22, cex=0.6, xpd = TRUE)
dev.off()


















