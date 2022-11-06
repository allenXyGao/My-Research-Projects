rm(list=ls())
library(survival)
library(readxl)
library(dplyr)
library(tidyr)


model = function(cov, geno) {
  res = coxph(Surv(days_os, PatientStatus_Description) ~ DxAge + BMI2 + as.factor(ER) + as.factor(PR) + as.factor(HER2) + as.factor(Hormonal.Therapy) +
                as.factor(Radiation.Therapy) + as.factor(Grade_Description) + as.factor(StageNEW) + as.factor(Surgery.Type) + geno, cbind(cov, geno=geno)) 
  print(res)
  return(summary(res)$coefficients["geno",])	
}

model.removed_HER2_Grade = function(cov, geno) {
  res = coxph(Surv(days_os, PatientStatus_Description) ~ DxAge + BMI2 + as.factor(ER) + as.factor(PR) + as.factor(Hormonal.Therapy) +
                as.factor(Radiation.Therapy) +  as.factor(StageNEW) + as.factor(Surgery.Type) + geno, cbind(cov, geno=geno)) 
  print(res)
  return(summary(res)$coefficients["geno",])	
} 

# no HER2 as covariate
model2 = function(cov, geno) {
  res = coxph(Surv(days_os, PatientStatus_Description) ~ DxAge + BMI2 + as.factor(ER) + as.factor(PR) + as.factor(Hormonal.Therapy) +
                as.factor(Radiation.Therapy) + as.factor(Grade_Description) + as.factor(StageNEW) + as.factor(Surgery.Type) + geno, cbind(cov, geno=geno))
  return(summary(res)$coefficients["geno",])
}

pheno = read.csv("/projects/rpci/qzhu/eschille/BreastGWAS/No_Y_Chrom/maf1.noHM/DBBR/Input/DBBR.pheno.csv", na.strings = c("NA", ''))
pheno = within(pheno, ER <- relevel(as.factor(ER), ref = "2"))
pheno = within(pheno, PR <- relevel(as.factor(PR), ref = "2"))
pheno = within(pheno, HER2 <- relevel(as.factor(HER2), ref = "2"))
pheno$DxAge = as.numeric(as.character(pheno$DxAge))
#pheno = within(pheno, Surgery.Type.clean <- relevel(as.factor(Surgery.Type.clean), ref = "None"))
pheno = within(pheno, Surgery.Type <- relevel(as.factor(Surgery.Type), ref = "None"))

pheno$rs11855431 = as.numeric(as.character(pheno$rs11855431))
pheno$rs6494889 = as.numeric(as.character(pheno$rs6494889))
pheno$rs28607477 = as.numeric(as.character(pheno$rs28607477))
pheno$rs720251 = as.numeric(as.character(pheno$rs720251))

#-------------------------------------------------------------------------------------------------------#
# step 1: merge 4 drugs with the oheno file
filepath.docetaxel <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/output_docetaxel/joined.docetaxel.ID_Date.xls"
filepath.doxorubicin <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/output_doxorubicin/out.doxorubicin.full_join.ID_Date.xls"
filepath.paclitaxol <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/output_taxol/updated.Taxol.ID_Date.only.xls"
filepath.trastuzumab <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/output_Trastuzumab/Trastuzumab.ID_Date.xls"

dat.docetaxel <- read_excel(path = filepath.docetaxel)
dat.doxorubicin <- read_excel(path = filepath.doxorubicin)
dat.paclitaxol <- read_excel(path = filepath.paclitaxol)
dat.trastuzumab <- read_excel(path = filepath.trastuzumab)

pheno.merged.drugs <- full_join(pheno, dat.docetaxel, by="SubjectID")
colnames(pheno.merged.drugs)[length(colnames(pheno.merged.drugs))]<- "date_docetaxel"
pheno.merged.drugs <- full_join(pheno.merged.drugs, dat.doxorubicin, by="SubjectID")
colnames(pheno.merged.drugs)[length(colnames(pheno.merged.drugs))] <- "date_doxorubicin"
pheno.merged.drugs <- full_join(pheno.merged.drugs, dat.paclitaxol, by=c("SubjectID"="SUBJECTID"))
colnames(pheno.merged.drugs)[length(colnames(pheno.merged.drugs))] <- "date_paclitaxol"
pheno.merged.drugs <- full_join(pheno.merged.drugs, dat.trastuzumab, by="SubjectID")
colnames(pheno.merged.drugs)[length(colnames(pheno.merged.drugs))] <- "date_trastuzumab"


#-------------------------------------------------------------------------------------------------------#

# step 2: create a new variable "date_par4_chemotherapy"

pheno.merged.drugs$date_par4_chemotherapy <- NA
for (i in 1:dim(pheno.merged.drugs)[1]) {
  date.doxorubicin <- pheno.merged.drugs$date_doxorubicin[i]
  date.trastuzumab <- pheno.merged.drugs$date_trastuzumab[i]
  if (!is.na(date.doxorubicin) & !is.na(date.trastuzumab)) {
    pheno.merged.drugs$date_par4_chemotherapy[i] <- min(date.doxorubicin, date.trastuzumab)
  } else if (!is.na(date.doxorubicin) & is.na(date.trastuzumab)) {
    pheno.merged.drugs$date_par4_chemotherapy[i] <- date.doxorubicin
  } else if (is.na(date.doxorubicin) & !is.na(date.trastuzumab))
    pheno.merged.drugs$date_par4_chemotherapy[i] <- date.trastuzumab
}


sum(is.na(pheno.merged.drugs$date_par4_chemotherapy)) # 3 removed, 467 left
#which(is.na(pheno.merged.drugs$date_par4_chemotherapy))

# please merge the dates you have for the 4 drugs as well as the 
# /projects/rpci/qzhu/eschille/BreastGWAS/No_Y_Chrom/maf1.noHM/DBBR/Input/DBBR.pheno.csv  file with this excel file ("Query1" sheet). 
# Please use "SubjectID" for merging across files 

filepath.query1 <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/inputs/select_patients/Yao Available Blood 11.1.21.xlsx"
query1.dat <- read_excel(path = filepath.query1, sheet = "Query1")
dim(query1.dat)
# 467   9
length(unique(query1.dat$SubjectID))
# 456
out.merged.drugs_pheno_query1 <- full_join(pheno.merged.drugs, query1.dat, by="SubjectID", suffix=c("_pheno_drugs","_query1"))
write.csv(out.merged.drugs_pheno_query1,
          file = "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/out.merged.drugs_pheno_query1.csv")




# received anthracycline or HER2 target therapy
newdat = pheno.merged.drugs[!is.na(pheno.merged.drugs$date_par4_chemotherapy),]
dim(newdat)

# write.csv(newdat, 
#           file = "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/test_Surv_associations/reveived_par4_chemotherapy.csv") 


#-------------------------------------------------------------------------------------------------------#

# step 3: 2*2 table of taxel vs docetaxel
newdat$flag_paclitaxol <- ifelse(!is.na(newdat$date_paclitaxol), "nonNA", "isNA")
newdat$flag_docetaxel <- ifelse(!is.na(newdat$date_docetaxel), "nonNA", "isNA")
table(newdat$flag_docetaxel, newdat$flag_paclitaxol)

#                     paclitaxol
#                     isNA   nonNA  Total
#             isNA     34     366    400
# docetaxel   nonNA    49     18     67
#             Total    83     384    467


#-------------------------------------------------------------------------------------------------------#

# step 4: test association of the 4 snps with OS in patients received par4_chemotherapy
cov = newdat[,c("days_os", "PatientStatus_Description", "DxAge", "BMI2", "ER", "PR", "HER2", "Hormonal.Therapy", "Radiation.Therapy", "Grade_Description", "StageNEW", "Surgery.Type")]
geno = newdat[,c("rs6494889", "rs28607477", "rs720251", "rs11855431")]
res = apply(geno, 2, function(x){model(cov=cov, geno=x)})
write.table(t(res), file="/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/test_Surv_associations/DBBR.task4.test", append=F, quote=F, sep="\t", col.names=T, row.names=T)	

# "rs6494889"
# total sample size=466 , n=411, number of events=67
# (55 observations deleted due to missingness)

# "rs28607477"
# total sample size=466, n= 414, number of events= 69 
# (52 observations deleted due to missingness)

# "rs720251"
# total sample size=466, n= 411, number of events= 68 
# (55 observations deleted due to missingness)

# "rs11855431"
# total sample size=466, n= 411, number of events= 68 
# (55 observations deleted due to missingness)

out <- c()
for (g in 1:4) {
  tmp <- cbind(cov, geno=geno[,g])
  out <- rbind(out, colSums(is.na(tmp)))
}
rownames(out) <- c("rs6494889", "rs28607477", "rs720251", "rs11855431")
write.table(out, file="/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/test_Surv_associations/DBBR.task4.test.missingness", append=F, quote=F, sep="\t", col.names=T, row.names=T)	
write.csv(out, file = "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/test_Surv_associations/DBBR.task4.test.missingness.csv") 

# removed 2 columns with missing: HER2 and Grade_Description
cov = newdat[,c("days_os", "PatientStatus_Description", "DxAge", "BMI2", "ER", "PR", "Hormonal.Therapy", "Radiation.Therapy",  "StageNEW", "Surgery.Type")]
geno = newdat[,c("rs6494889", "rs28607477", "rs720251", "rs11855431")]
res = apply(geno, 2, function(x){model.removed_HER2_Grade(cov=cov, geno=x)})
write.table(t(res), file="/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/test_Surv_associations/DBBR.task4.removed_HER2_Grade.test", append=F, quote=F, sep="\t", col.names=T, row.names=T)	



# "rs6494889"
# total sample size=466, n= 451, number of events= 77 
# (15 observations deleted due to missingness)

# "rs28607477"
# n= 454, number of events= 79 
# (12 observations deleted due to missingness)

# "rs720251"
# n= 451, number of events= 78 
# (15 observations deleted due to missingness)

# "rs11855431"
# n= 451, number of events= 78 
# (15 observations deleted due to missingness)


#-------------------------------------------------------------------------------------------------------#

# task 5: two groups: 1. received paclitaxel before par4_chemotherapy
#                     2. otherwise

newdat$paclitaxel_earlier_par4chemotherapy <- NA
for (i in 1:dim(newdat)[1]) {
  date.paclitaxel <- newdat$date_paclitaxol[i]
  date.par4chemotherapy <- newdat$date_par4_chemotherapy[i]
  if (!is.na(date.par4chemotherapy) & !is.na(date.paclitaxel)) {
    if (date.paclitaxel < date.par4chemotherapy) {
      newdat$paclitaxel_earlier_par4chemotherapy[i] <- "YES"
    } else if (date.paclitaxel == date.par4chemotherapy) {
      newdat$paclitaxel_earlier_par4chemotherapy[i] <- "EQUAL"
    } else  {
      newdat$paclitaxel_earlier_par4chemotherapy[i] <- "Later"
    }
  } else if (!is.na(date.par4chemotherapy) | !is.na(date.paclitaxel)) {
    newdat$paclitaxel_earlier_par4chemotherapy[i] <- "hasNA"
  }
}


dat.paclitaxel_earlier_par4chemotherapy <- newdat[newdat$paclitaxel_earlier_par4chemotherapy == "YES",]

write.csv(dat.paclitaxel_earlier_par4chemotherapy, 
          file = "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/test_Surv_associations/paclitaxel_earlier_par4chemotherapy.csv") 

#-------------------------------------------------------------------------------------------------------#
# task 5.1
# could you put the 8 PT who received paclitaxel before par4-dependent chemo therapy and 
# the 34 pts who didn't receive paclitaxel or docetaxel in the same group and 
# test association between the SNP genotypes and OS while controlling for all covariates? 
# please do the same for patients not in this group.

group.index <- c(which(newdat$paclitaxel_earlier_par4chemotherapy == "YES"), 
           which( (newdat$flag_paclitaxol == "isNA") & (newdat$flag_docetaxel == "isNA") ))
# PT in the group
cov = newdat[group.index,c("days_os", "PatientStatus_Description", "DxAge", "BMI2", "ER", "PR", "HER2", "Hormonal.Therapy", "Radiation.Therapy", "Grade_Description", "StageNEW", "Surgery.Type")]
geno = newdat[group.index,c("rs6494889", "rs28607477", "rs720251", "rs11855431")]

# "rs6494889"
tmp <- cbind(cov, geno=geno[,1])
tmp <- tmp[complete.cases(tmp), ]
apply(tmp, 2, table)
tmp$Grade_Description[tmp$Grade_Description == 1] <- 2
res = coxph(Surv(days_os, PatientStatus_Description) ~ DxAge + BMI2 + as.factor(ER) + as.factor(PR) + as.factor(HER2) + as.factor(Hormonal.Therapy) +
              as.factor(Radiation.Therapy) + as.factor(Grade_Description) + as.factor(StageNEW) + as.factor(Surgery.Type) + geno, tmp) 
print(res)
summary(res)$coefficients["geno",]
# coef   exp(coef)    se(coef)           z    Pr(>|z|) 
# -3.59270371  0.02752381  2.57050947 -1.39766212  0.16221456 
# sample size= 35, number of events= 9 

# "rs28607477"
tmp <- cbind(cov, geno=geno[,2])
tmp <- tmp[complete.cases(tmp), ]
apply(tmp, 2, table)
tmp$Grade_Description[tmp$Grade_Description == 1] <- 2
res = coxph(Surv(days_os, PatientStatus_Description) ~ DxAge + BMI2 + as.factor(ER) + as.factor(PR) + as.factor(HER2) + as.factor(Hormonal.Therapy) +
              as.factor(Radiation.Therapy) + as.factor(Grade_Description) + as.factor(StageNEW) + as.factor(Surgery.Type) + geno, tmp) 
print(res)
summary(res)$coefficients["geno",]
# coef   exp(coef)    se(coef)           z    Pr(>|z|) 
# -3.40859445  0.03308767  2.50825285 -1.35895169  0.17416190 
# sample size= 35, number of events= 9 

# "rs720251"
tmp <- cbind(cov, geno=geno[,3])
tmp <- tmp[complete.cases(tmp), ]
apply(tmp, 2, table)
tmp$Grade_Description[tmp$Grade_Description == 1] <- 2
res = coxph(Surv(days_os, PatientStatus_Description) ~ DxAge + BMI2 + as.factor(ER) + as.factor(PR) + as.factor(HER2) + as.factor(Hormonal.Therapy) +
              as.factor(Radiation.Therapy) + as.factor(Grade_Description) + as.factor(StageNEW) + as.factor(Surgery.Type) + geno, tmp) 
print(res)
summary(res)$coefficients["geno",]
# coef   exp(coef)    se(coef)           z    Pr(>|z|) 
# -3.40859445  0.03308767  2.50825285 -1.35895169  0.17416190 
# sample size= 35, number of events= 9

# "rs11855431"
tmp <- cbind(cov, geno=geno[,4])
tmp <- tmp[complete.cases(tmp), ]
apply(tmp, 2, table)
tmp$Grade_Description[tmp$Grade_Description == 1] <- 2
res = coxph(Surv(days_os, PatientStatus_Description) ~ DxAge + BMI2 + as.factor(ER) + as.factor(PR) + as.factor(HER2) + as.factor(Hormonal.Therapy) +
              as.factor(Radiation.Therapy) + as.factor(Grade_Description) + as.factor(StageNEW) + as.factor(Surgery.Type) + geno, tmp) 
print(res)
summary(res)$coefficients["geno",]
# coef   exp(coef)    se(coef)           z    Pr(>|z|) 
# -3.59270371  0.02752381  2.57050947 -1.39766212  0.16221456 
# sample size= 35, number of events= 9

cov = newdat[group.index, c("days_os", "PatientStatus_Description", "DxAge", "BMI2", "ER", "PR", 
                            "Hormonal.Therapy", "Radiation.Therapy", "StageNEW", "Surgery.Type")]
geno = newdat[group.index, c("rs6494889", "rs28607477", "rs720251", "rs11855431")]
res = apply(geno, 2, function(x){model.removed_HER2_Grade(cov=cov, geno=x)})
write.table(t(res), file="/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/test_Surv_associations/DBBR.8.Before_34.no_paclitaxel.docetaxel_removed_HER2.GRADE", append=F, quote=F, sep="\t", col.names=T, row.names=T)	

# "rs6494889"
# total sample size=42, n= 40, number of events= 13 
# (2 observations deleted due to missingness)

# "rs28607477"
# total sample size=42, n= 40, number of events= 13 
# (2 observations deleted due to missingness)

# "rs720251"
# total sample size=42, n= 40, number of events= 13 
# (2 observations deleted due to missingness)

# "rs11855431"
# total sample size=42, n= 40, number of events= 13 
# (2 observations deleted due to missingness)

# PT not in the group

cov = newdat[-group.index,c("days_os", "PatientStatus_Description", "DxAge", "BMI2", "ER", "PR", "HER2", "Hormonal.Therapy", "Radiation.Therapy", "Grade_Description", "StageNEW", "Surgery.Type")]
geno = newdat[-group.index,c("rs6494889", "rs28607477", "rs720251", "rs11855431")]
res = apply(geno, 2, function(x){model(cov=cov, geno=x)})
write.table(t(res), file="/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/test_Surv_associations/DBBR.NOT_8.Before_34.no_paclitaxel.docetaxel", append=F, quote=F, sep="\t", col.names=T, row.names=T)	
# "rs6494889"
# total sample size=424, n= 376, number of events= 58 
#  (48 observations deleted due to missingness)

# "rs28607477"
# total sample size=424,  n= 379, number of events= 60 
#  (45 observations deleted due to missingness)

# "rs720251"
# total sample size=424, n= 376, number of events= 59 
#  (48 observations deleted due to missingness)

# "rs11855431"
# total sample size=424, n= 376, number of events= 59 
#   (48 observations deleted due to missingness)


cov = newdat[-group.index, c("days_os", "PatientStatus_Description", "DxAge", "BMI2", "ER", "PR", 
                            "Hormonal.Therapy", "Radiation.Therapy", "StageNEW", "Surgery.Type")]
geno = newdat[-group.index, c("rs6494889", "rs28607477", "rs720251", "rs11855431")]
res = apply(geno, 2, function(x){model.removed_HER2_Grade(cov=cov, geno=x)})
write.table(t(res), file="/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/test_Surv_associations/DBBR.NOT_8.Before_34.no_paclitaxel.docetaxel_removed_HER2.GRADE", append=F, quote=F, sep="\t", col.names=T, row.names=T)	

# "rs6494889"
# total sample size=424, n= 411, number of events= 64 
# (13 observations deleted due to missingness)

# "rs28607477"
# total sample size=424,  n= 414, number of events= 66 
# (10 observations deleted due to missingness)

# "rs720251"
# total sample size=424, n= 411, number of events= 65 
# (13 observations deleted due to missingness)

# "rs11855431"
# total sample size=424, n= 411, number of events= 65 
# (13 observations deleted due to missingness)


#-------------------------------------------------------------------------------------------------------#
# could you put the 8 PT who received paclitaxel before par4-dependent chemo therapy and 
# the 82 pts who didn't receive paclitaxel in the same group and test the associations? 
# pls do the same for the patients not in this group

group.index <- c(which(newdat$paclitaxel_earlier_par4chemotherapy == "YES"), 
                 which(newdat$flag_paclitaxol == "isNA" ))

# PT in the group
cov = newdat[group.index,c("days_os", "PatientStatus_Description", "DxAge", "BMI2", "ER", "PR", "HER2", "Hormonal.Therapy", "Radiation.Therapy", "Grade_Description", "StageNEW", "Surgery.Type")]
geno = newdat[group.index,c("rs6494889", "rs28607477", "rs720251", "rs11855431")]
res = apply(geno, 2, function(x){model(cov=cov, geno=x)})
write.table(t(res), file="/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/test_Surv_associations/DBBR.8.Before_82.no_paclitaxel", append=F, quote=F, sep="\t", col.names=T, row.names=T)	

# "rs6494889"
# total sample size=90, n= 81, number of events= 15 
#  (9 observations deleted due to missingness)

# "rs28607477"
# total sample size=90, n= 82, number of events= 16 
# (8 observations deleted due to missingness)

# "rs720251"
# total sample size=90, n= 81, number of events= 15 
# (9 observations deleted due to missingness)

# "rs11855431"
# total sample size=90, n= 81, number of events= 15 
# (9 observations deleted due to missingness)


cov = newdat[group.index, c("days_os", "PatientStatus_Description", "DxAge", "BMI2", "ER", "PR", 
                            "Hormonal.Therapy", "Radiation.Therapy", "StageNEW", "Surgery.Type")]
geno = newdat[group.index, c("rs6494889", "rs28607477", "rs720251", "rs11855431")]
res = apply(geno, 2, function(x){model.removed_HER2_Grade(cov=cov, geno=x)})
write.table(t(res), file="/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/test_Surv_associations/DBBR.8.Before_82.no_paclitaxel_removed_HER2.GRADE", append=F, quote=F, sep="\t", col.names=T, row.names=T)	

# "rs6494889"
# total sample size=90, n= 86, number of events= 19 
# (4 observations deleted due to missingness)

# "rs28607477"
# total sample size=90, n= 87, number of events= 20 
# (3 observations deleted due to missingness)

# "rs720251"
# total sample size=90, n= 86, number of events= 19 
# (4 observations deleted due to missingness)

# "rs11855431"
# total sample size=90, n= 86, number of events= 19 
# (4 observations deleted due to missingness)


# PT not in the group
cov = newdat[-group.index,c("days_os", "PatientStatus_Description", "DxAge", "BMI2", "ER", "PR", "HER2", "Hormonal.Therapy", "Radiation.Therapy", "Grade_Description", "StageNEW", "Surgery.Type")]
geno = newdat[-group.index,c("rs6494889", "rs28607477", "rs720251", "rs11855431")]
res = apply(geno, 2, function(x){model(cov=cov, geno=x)})
write.table(t(res), file="/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/test_Surv_associations/DBBR.NOT_8.Before_82.no_paclitaxel", append=F, quote=F, sep="\t", col.names=T, row.names=T)	

# "rs6494889"
# total sample size=376, n= 330, number of events= 52 
# (46 observations deleted due to missingness)

# "rs28607477"
# total sample size=376, n= 332, number of events= 53 
# (44 observations deleted due to missingness)

# "rs720251"
# total sample size=376, n= 330, number of events= 53 
# (46 observations deleted due to missingness)

# "rs11855431"
# total sample size=376, n= 330, number of events= 53 
# (46 observations deleted due to missingness)


cov = newdat[-group.index, c("days_os", "PatientStatus_Description", "DxAge", "BMI2", "ER", "PR", 
                            "Hormonal.Therapy", "Radiation.Therapy", "StageNEW", "Surgery.Type")]
geno = newdat[-group.index, c("rs6494889", "rs28607477", "rs720251", "rs11855431")]
res = apply(geno, 2, function(x){model.removed_HER2_Grade(cov=cov, geno=x)})
write.table(t(res), file="/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/test_Surv_associations/DBBR.NOT_8.Before_82.no_paclitaxel_removed_HER2.GRADE", append=F, quote=F, sep="\t", col.names=T, row.names=T)	

# "rs6494889"
# total sample size=376, n= 365, number of events= 58 
# (11 observations deleted due to missingness)

# "rs28607477"
# total sample size=376, n= 367, number of events= 59 
# (9 observations deleted due to missingness)

# "rs720251"
# total sample size=376, n= 365, number of events= 59 
# (11 observations deleted due to missingness)

# "rs11855431"
# total sample size=376, n= 365, number of events= 59 
# (11 observations deleted due to missingness)





#-------------------------------------------------------------------------------------------------------#

# could you first test the association of the 4 SNPs with OS in patients with Paclitaxel and patients without paclitaxel?


table(newdat$flag_paclitaxol)
# isNA nonNA 
#  82   384 
newdat.with_paclitaxol <- newdat %>% filter(flag_paclitaxol == "nonNA")
newdat.without_paclitaxol <- newdat %>% filter(flag_paclitaxol == "isNA")

# with paclitaxol
cov = newdat.with_paclitaxol[,c("days_os", "PatientStatus_Description", "DxAge", "BMI2", "ER", "PR", "HER2", "Hormonal.Therapy", "Radiation.Therapy", "Grade_Description", "StageNEW", "Surgery.Type")]
geno = newdat.with_paclitaxol[,c("rs6494889", "rs28607477", "rs720251", "rs11855431")]
res = apply(geno, 2, function(x){model(cov=cov, geno=x)})
write.table(t(res), file="/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/test_Surv_associations/DBBR.with.paclitaxol", append=F, quote=F, sep="\t", col.names=T, row.names=T)	


# "rs6494889"
# total sample size=384, n= 336, number of events= 55 
# (48 observations deleted due to missingness)

# "rs28607477"
# total sample size=384, n= 338, number of events= 56 
# (46 observations deleted due to missingness)

# "rs720251"
# total sample size=384, n= 336, number of events= 56 
# (48 observations deleted due to missingness)

# "rs11855431"
# total sample size=384, n= 336, number of events= 56 
# (48 observations deleted due to missingness)


# removed her2 + grade
cov = newdat.with_paclitaxol[,c("days_os", "PatientStatus_Description", "DxAge", "BMI2", "ER", "PR", "Hormonal.Therapy", "Radiation.Therapy",  "StageNEW", "Surgery.Type")]
geno = newdat.with_paclitaxol[,c("rs6494889", "rs28607477", "rs720251", "rs11855431")]
res = apply(geno, 2, function(x){model.removed_HER2_Grade(cov=cov, geno=x)})
write.table(t(res), file="/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/test_Surv_associations/DBBR.with.paclitaxol.removed_HER2_Grade.test", append=F, quote=F, sep="\t", col.names=T, row.names=T)	

# "rs6494889"
# total sample size=384, n= 373, number of events= 63 
# (11 observations deleted due to missingness)

# "rs28607477"
# total sample size=384,n= 375, number of events= 64 
# (9 observations deleted due to missingness)

# "rs720251"
# total sample size=384, n= 373, number of events= 64 
# (11 observations deleted due to missingness)

# "rs11855431"
# total sample size=384, n= 373, number of events= 64 
# (11 observations deleted due to missingness)





# without paclitaxol
cov = newdat.without_paclitaxol[,c("days_os", "PatientStatus_Description", "DxAge", "BMI2", "ER", "PR", "HER2", "Hormonal.Therapy", "Radiation.Therapy", "Grade_Description", "StageNEW", "Surgery.Type")]
geno = newdat.without_paclitaxol[,c("rs6494889", "rs28607477", "rs720251", "rs11855431")]
res = apply(geno, 2, function(x){model(cov=cov, geno=x)})
write.table(t(res), file="/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/test_Surv_associations/DBBR.without.paclitaxol", append=F, quote=F, sep="\t", col.names=T, row.names=T)	

# "rs6494889"
# total sample size=82, n= 75, number of events= 12 
# (7 observations deleted due to missingness)

# "rs28607477"
# total sample size=82, n= 76, number of events= 13 
# (6 observations deleted due to missingness)

# "rs720251"
# total sample size=82, n= 75, number of events= 12 
# (7 observations deleted due to missingness)

# "rs11855431"
# total sample size=82,n= 75, number of events= 12 
# (7 observations deleted due to missingness)

# removed HER2 + Grade
cov = newdat.without_paclitaxol[,c("days_os", "PatientStatus_Description", "DxAge", "BMI2", "ER", "PR", "Hormonal.Therapy", "Radiation.Therapy", "StageNEW", "Surgery.Type")]
geno = newdat.without_paclitaxol[,c("rs6494889", "rs28607477", "rs720251", "rs11855431")]
res = apply(geno, 2, function(x){model.removed_HER2_Grade(cov=cov, geno=x)})
write.table(t(res), file="/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/test_Surv_associations/DBBR.without.paclitaxol.removed_HER2_Grade.test", append=F, quote=F, sep="\t", col.names=T, row.names=T)	

# "rs6494889"
# total sample size=82,n= 78, number of events= 14 
# (4 observations deleted due to missingness)

# "rs28607477"
# total sample size=82,n= 79, number of events= 15 
# (3 observations deleted due to missingness)

# "rs720251"
# total sample size=82, n= 78, number of events= 14 
#(4 observations deleted due to missingness)

# "rs11855431"
# total sample size=82,n= 78, number of events= 14 
#(4 observations deleted due to missingness)



#-------------------------------------------------------------------------------------------------------
# please do the same for with docetaxel and without docetaxel

table(newdat$flag_docetaxel)
# isNA nonNA 
# 400    66
newdat.with_docetaxel <- newdat %>% filter(flag_docetaxel == "nonNA")
newdat.without_docetaxel <- newdat %>% filter(flag_docetaxel == "isNA")

# with docetaxel
cov = newdat.with_docetaxel[,c("days_os", "PatientStatus_Description", "DxAge", "BMI2", "ER", "PR", "HER2", "Hormonal.Therapy", "Radiation.Therapy", "Grade_Description", "StageNEW", "Surgery.Type")]
geno = newdat.with_docetaxel[,c("rs6494889", "rs28607477", "rs720251", "rs11855431")]
res = apply(geno, 2, function(x){model(cov=cov, geno=x)})
write.table(t(res), file="/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/test_Surv_associations/DBBR.with.docetaxel", append=F, quote=F, sep="\t", col.names=T, row.names=T)	

# "rs6494889"
# total sample size=66, n= 61, number of events= 10 
# (5 observations deleted due to missingness)

# "rs28607477"
# total sample size=66,n= 62, number of events= 11 
# (4 observations deleted due to missingness)

# "rs720251"
# total sample size=66,n= 61, number of events= 10 
# (5 observations deleted due to missingness)

# "rs11855431"
# total sample size=66,n= 61, number of events= 10 
# (5 observations deleted due to missingness)

# removed HER2+ grade
cov = newdat.with_docetaxel[,c("days_os", "PatientStatus_Description", "DxAge", "BMI2", "ER", "PR",  "Hormonal.Therapy", "Radiation.Therapy",  "StageNEW", "Surgery.Type")]
geno = newdat.with_docetaxel[,c("rs6494889", "rs28607477", "rs720251", "rs11855431")]
res = apply(geno, 2, function(x){model.removed_HER2_Grade(cov=cov, geno=x)})
write.table(t(res), file="/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/test_Surv_associations/DBBR.with.docetaxel.removed_HER2_Grade", append=F, quote=F, sep="\t", col.names=T, row.names=T)	

# "rs6494889"
# total sample size=66, n= 62, number of events= 10 
# (4 observations deleted due to missingness)

# "rs28607477"
# total sample size=66,n= 63, number of events= 11 
# (3 observations deleted due to missingness)

# "rs720251"
# total sample size=66,n= 62, number of events= 10 
# (4 observations deleted due to missingness)

# "rs11855431"
# total sample size=66,n= 62, number of events= 10 
# (4 observations deleted due to missingness)



# without docetaxel
cov = newdat.without_docetaxel[,c("days_os", "PatientStatus_Description", "DxAge", "BMI2", "ER", "PR", "HER2", "Hormonal.Therapy", "Radiation.Therapy", "Grade_Description", "StageNEW", "Surgery.Type")]
geno = newdat.without_docetaxel[,c("rs6494889", "rs28607477", "rs720251", "rs11855431")]
res = apply(geno, 2, function(x){model(cov=cov, geno=x)})
write.table(t(res), file="/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/test_Surv_associations/DBBR.without.docetaxel", append=F, quote=F, sep="\t", col.names=T, row.names=T)	

# "rs6494889"
# total sample size=400, n= 350, number of events= 57 
# (50 observations deleted due to missingness)
# 
# "rs28607477"
# total sample size=400, n= 352, number of events= 58 
# (48 observations deleted due to missingness)
# 
# "rs720251"
# total sample size=400, n= 350, number of events= 58 
# (50 observations deleted due to missingness)
# 
# "rs11855431"
# total sample size=400, n= 350, number of events= 58 
# (50 observations deleted due to missingness)

# removed HER2+grade
cov = newdat.without_docetaxel[,c("days_os", "PatientStatus_Description", "DxAge", "BMI2", "ER", "PR", "Hormonal.Therapy", "Radiation.Therapy", "StageNEW", "Surgery.Type")]
geno = newdat.without_docetaxel[,c("rs6494889", "rs28607477", "rs720251", "rs11855431")]
res = apply(geno, 2, function(x){model.removed_HER2_Grade(cov=cov, geno=x)})
write.table(t(res), file="/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/test_Surv_associations/DBBR.without.docetaxel.removed_HER2_Grade", append=F, quote=F, sep="\t", col.names=T, row.names=T)	


# "rs6494889"
# total sample size=400, n= 389, number of events= 67 
# (11 observations deleted due to missingness)
# 
# "rs28607477"
# total sample size=400, n= 391, number of events= 68 
# (9 observations deleted due to missingness)
# 
# "rs720251"
# total sample size=400, n= 389, number of events= 68 
# (11 observations deleted due to missingness)
# 
# "rs11855431"
# total sample size=400, n= 389, number of events= 68 
# (11 observations deleted due to missingness)


