# create phenotype data files
rm(list = ls())
library(dplyr)
library(tidyr)

pop <- c("ALLcasecontrol_cohort1", "ALLcasecontrol_cohort2",
         "AMLMDScasecontrol_cohort1", "AMLMDScasecontrol_cohort2")


i=3



#### PC
filepath <- paste("/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/EigenStrat_Results/",
                  pop[i],"/", sep="")
PCs <- read.table(file=paste(filepath, "LDprune_autosome-pca2.eth.out", sep=""), header = T)
head(PCs)
PC.names <- c()
for (i in 1:20) {
  PC.names <- c(PC.names, paste0("PC", i))
}
colnames(PCs)[3:dim(PCs)[2]] <- PC.names
head(PCs)

# extract significant PCs for this specific file
twout <- read.table(file=paste(filepath, "LDprune_autosome-pca/LDprune_autosome-pca2.twout", sep=""))
sig.n <- sum(twout[1:20,5]<=0.05)  ## use top 20 to calculate, sig.n<20
print(sig.n)

PCs.sig.cut <- PCs[, c(3:(3+sig.n-1))]
head(PCs.sig.cut)

# identify FID/IID
out <- strsplit(as.character(PCs[, 1]),':')
FID_IID.remaining <- do.call(rbind, out)
PCs.sig.FID.IID <- cbind(FID_IID.remaining, PCs.sig.cut)
colnames(PCs.sig.FID.IID)[1:2] <- c("FID", "short_IID")
head(PCs.sig.FID.IID)
#----------------------------------------------------------------------------------------------


#### Sample ID + age + sex
filepath <- "/projects/rpci/qzhu/lsuchest/lsuchest/DBMT_PhenoData/"
#meta_phenotype <- read.table(file=paste(filepath, "DBMT_PhenoData_EA_long_allVar_20190223.txt", sep=""), fill = T)
meta_phenotype <- read.table(file=paste(filepath, "DBMT_PhenoData_EA_long_allVar_20190223.txt", sep=""), sep="\t", header = T)
meta_phenotype.less <- meta_phenotype[, c("Exome.FID", "Exome.IID", "Exome.recip.pop", "Exome.donor.pop", 
                                     "age", "dnrage", "sex", "dnrsex")]
head(meta_phenotype.less)
dim(meta_phenotype.less)
#----------------------------------------------------------------------------------------------


#### case/control
ALL.case.cohort1 <- read.table(file="/projects/rpci/qzhu/eschille/transplant.GRANT/data/Recipients/EAcohort1/ALL.samples.txt")
# cases=2, control=1
ALL.case.cohort1$phenotype <- 2
head(ALL.case.cohort1)
control.cohort1 <- read.table(file="/projects/rpci/qzhu/qzhu/mywork/LaraBMTexome/input/BMT.snpQC.final.donor.EA.cohort1.fam")[, 1:2]
control.cohort1$phenotype <- 1
head(control.cohort1)
ALL.casecontrol.cohort1 <- rbind(control.cohort1, ALL.case.cohort1)
colnames(ALL.casecontrol.cohort1)[1:2] <- c("FID", "IID")
head(ALL.casecontrol.cohort1)
#----------------------------------------------------------------------------------------------




### merge data sets
# merge PCs.sig.FID.IID and ALL.casecontrol.cohort1 -> res1
out2 <- strsplit(ALL.casecontrol.cohort1$IID, "([A-Z]-)", perl = TRUE)
ALL.casecontrol.cohort1$short_IID <- do.call(rbind, out2)[,2]
ALL.casecontrol.cohort1$FID <- as.character(ALL.casecontrol.cohort1$FID)
ALL.casecontrol.cohort1$short_IID <- as.character(ALL.casecontrol.cohort1$short_IID)
PCs.sig.FID.IID$FID <- as.character(PCs.sig.FID.IID$FID)
PCs.sig.FID.IID$short_IID <- as.character(PCs.sig.FID.IID$short_IID)
res1 <- ALL.casecontrol.cohort1 %>% inner_join(PCs.sig.FID.IID, by=c("FID", "short_IID"))

### Problems
# merge res1 and meta_phenotype.less
#out3 <- strsplit(meta_phenotype.less$Exome.IID, "([A-Z]-)", perl = TRUE)
#meta_phenotype.less$short_IID <- do.call(rbind, out3)[,2]
#meta_phenotype.less$short_IID <- as.character(meta_phenotype.less$short_IID)

meta_phenotype.less$Exome.FID <- as.character(meta_phenotype.less$Exome.FID)
meta_phenotype.less$Exome.IID <- as.character(meta_phenotype.less$Exome.IID)

colnames(meta_phenotype.less)[colnames(meta_phenotype.less) == "Exome.FID"] <- "FID"
colnames(meta_phenotype.less)[colnames(meta_phenotype.less) == "Exome.IID"] <- "IID"

### here I used inner_join
# --------------- Problems -----------------
res2 <- res1 %>% inner_join(meta_phenotype.less, by=c("FID", "IID"))

#res.try <- res1 %>% left_join(meta_phenotype.less, by=c("FID", "IID"))
# --------------- Problems -----------------
head(res2 )


#### after merge, select columns
# cases=2, control=1
res2$data.sex <- NA
res2$data.sex[res2$phenotype == 1] <-  res2[res2$phenotype == 1, "dnrsex"]
res2$data.sex[res2$phenotype == 2] <-  res2[res2$phenotype == 2, "sex"]
res2$data.age <- NA
res2$data.age[res2$phenotype == 1] <- res2[res2$phenotype == 1, "dnrage"]
res2$data.age[res2$phenotype == 2] <- res2[res2$phenotype == 2, "age"]

res.phenotype <- res2[, c("FID", "IID", "phenotype", "data.age", "data.sex", PC.names[1:sig.n])]
colnames(res.phenotype)[colnames(res.phenotype) == "data.age"] = "age"
colnames(res.phenotype)[colnames(res.phenotype) == "data.sex"] = "sex"
head(res.phenotype)

#-----------------------------------------------
#### cases
# ALL (cases in cohort 1): /projects/rpci/qzhu/eschille/transplant.GRANT/data/Recipients/EAcohort1/ALL.samples.txt

# AML.MDS(cases in cohort 1): /projects/rpci/qzhu/eschille/transplant.GRANT/data/Recipients/EAcohort1/AML.MDS.samples.txt

# ALL(cases in cohort 2): /projects/rpci/qzhu/eschille/transplant.GRANT/data/Recipients/EAcohort2/ALL.samples.txt

# AML.MDS(cases in cohort 2): /projects/rpci/qzhu/eschille/transplant.GRANT/data/Recipients/EAcohort2/AML.MDS.samples.txt


#### control

# (control in cohort 1): /projects/rpci/qzhu/qzhu/mywork/LaraBMTexome/input/BMT.snpQC.final.donor.EA.cohort1.fam

# (control in cohort 2): /projects/rpci/qzhu/qzhu/mywork/LaraBMTexome/input/BMT.snpQC.final.donor.EA.cohort2.fam





c1 <- dim(read.table(file="/projects/rpci/qzhu/eschille/transplant.GRANT/data/Recipients/EAcohort1/ALL.samples.txt"))[1]
c2 <- dim(read.table(file="/projects/rpci/qzhu/eschille/transplant.GRANT/data/Recipients/EAcohort1/AML.MDS.samples.txt"))[1]

c3 <- dim(read.table(file="/projects/rpci/qzhu/eschille/transplant.GRANT/data/Recipients/EAcohort2/ALL.samples.txt"))[1]
c4 <- dim(read.table(file="/projects/rpci/qzhu/eschille/transplant.GRANT/data/Recipients/EAcohort2/AML.MDS.samples.txt"))[1]

c5 <- dim(read.table(file="/projects/rpci/qzhu/qzhu/mywork/LaraBMTexome/input/BMT.snpQC.final.donor.EA.cohort1.fam"))[1]
c6 <- dim(read.table(file="/projects/rpci/qzhu/qzhu/mywork/LaraBMTexome/input/BMT.snpQC.final.donor.EA.cohort2.fam"))[1]
data.names <- c("ALL_case_C1", "AMLMDS_case_C1", "ALL_case_C2", "AMLMDS_case_C2", "Control_C1", "Control_C2")

ct <- c(c1,c2,c3,c4,c5,c6)
data.frame(rbind(data.names, ct))
total.ct <- c1+c2+c3+c4+c5+c6
total.ct


colSums((!is.na(meta_phenotype.less)))

