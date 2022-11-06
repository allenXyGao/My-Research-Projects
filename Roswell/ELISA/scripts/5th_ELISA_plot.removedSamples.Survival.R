rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)
# library(hrbrthemes)

filepath <- "/projects/rpci/qzhu/users-qzhu/grant/UACA_R01/ELISA/"
PT.genotype <- read.table(paste0(filepath, "5th_ELISA.select_samples.txt"), header = TRUE, sep = "\t")
PT.genotype.rm <- PT.genotype[- which(PT.genotype$CollectionID=="Cl-00044311" & PT.genotype$COLLECTIONBARCODE == "1002201"), ]

out.PT.genotype <- PT.genotype.rm %>% select(SubjectID, rs6494889)
dat.450_570 <- read.table(paste0(filepath, "5thELISA_450nm-570nm"), header = TRUE, colClasses = c("numeric", "numeric", "numeric"))
dat.450_540 <- read.table(paste0(filepath, "5thELISA_450nm-540nm"), header = TRUE, colClasses = c("numeric", "numeric", "numeric"))


ID.info <- readxl::read_excel("/projects/rpci/qzhu/users-qzhu/grant/UACA_R01/ELISA/DistSamplesList_nowaka_Sep-16-2022-11-25-05.xls") # collection_alias
ID.info <- ID.info %>% select(Subjectid, `Collection Alias`, Container) 
colnames(ID.info) <- c("Subjectid", "ID_alias", "Container")
ID.info <- unique( ID.info )
ID.info$ID_alias <- as.numeric(ID.info$ID_alias)

dat.450_540 <- left_join(dat.450_540, ID.info, by=c("ID"="ID_alias"))  
dat.450_570 <- left_join(dat.450_570, ID.info, by=c("ID"="ID_alias"))
length(unique(dat.450_570$Subjectid)) # 39 PTs
length(unique(dat.450_570$ID)) # 40 records


dat.450_540.genotype <- inner_join(dat.450_540, out.PT.genotype, by=c("Subjectid"="SubjectID"))
dat.450_570.genotype <- inner_join(dat.450_570, out.PT.genotype, by=c("Subjectid"="SubjectID"))
length(unique(dat.450_540.genotype$Subjectid))  # 38
length(unique(dat.450_570.genotype$Subjectid))  # 38 

dat.450_540.genotype$cleaved_par4 <- dat.450_540.genotype$original - 0.5*dat.450_540.genotype$X30K
dat.450_570.genotype$cleaved_par4 <- dat.450_570.genotype$original - 0.5*dat.450_570.genotype$X30K



# 11/4/2022 TASK:
# coxph(Surv(days_os, PatientStatus_Description) ~ DxAge + BMI2 + as.factor(ER) + as.factor(PR) + as.factor(Hormonal.Therapy) +
#         as.factor(Radiation.Therapy) + as.factor(Grade_Description) + as.factor(StageNEW) + as.factor(Surgery.Type) + Par-4 level 
#       for white straw and green straw seperately



head(PT.genotype.rm)
colSums(is.na(PT.genotype.rm))
pheno <- PT.genotype.rm

pheno = within(pheno, ER <- relevel(as.factor(ER), ref = "2"))
pheno = within(pheno, PR <- relevel(as.factor(PR), ref = "2"))
pheno = within(pheno, HER2 <- relevel(as.factor(HER2), ref = "2"))
pheno$DxAge = as.numeric(as.character(pheno$DxAge))
#pheno = within(pheno, Surgery.Type.clean <- relevel(as.factor(Surgery.Type.clean), ref = "None"))
pheno = within(pheno, Surgery.Type <- relevel(as.factor(Surgery.Type), ref = "Lumpectomy"))
pheno$rs6494889 = as.numeric(as.character(pheno$rs6494889))
pheno$datediff <- difftime(as.Date(pheno$COLLECTIONDATE, "%Y-%m-%d"), as.Date(pheno$DateDx, format = "%m/%d/%Y"), units = "days")


dat.450_540.merged <- dat.450_540.genotype %>% select(-rs6494889) 
dat.450_570.merged <- dat.450_570.genotype %>% select(-rs6494889)

pheno.450_540 <- inner_join(pheno, dat.450_540.merged, by=c("SubjectID"="Subjectid"))
pheno.450_570 <- inner_join(pheno, dat.450_570.merged, by=c("SubjectID"="Subjectid"))


library(survival)
# follow this script: /projects/rpci/qzhu/qzhu/mywork/PathwaysGWAS/scripts/DBBR.taxol.R

# 
# tmp_all <- dat%>% select(days_os, PatientStatus_Description, DxAge , BMI2 ,ER , PR,  Hormonal.Therapy,
#                                                                     Radiation.Therapy ,Grade_Description , StageNEW, Surgery.Type)
# tmp_green <- dat%>% filter(Container == "CBS Straw Green") %>% select(days_os, PatientStatus_Description, DxAge , BMI2 ,ER , PR,  Hormonal.Therapy,
#                                                                 Radiation.Therapy ,Grade_Description , StageNEW, Surgery.Type)
# tmp_white <- dat%>% filter(Container == "CBS Straw White") %>% select(days_os, PatientStatus_Description, DxAge , BMI2 ,ER , PR,  Hormonal.Therapy,
#                                                                       Radiation.Therapy ,Grade_Description , StageNEW, Surgery.Type)

dat = pheno.450_570
# white
res = coxph(Surv(days_os, PatientStatus_Description) ~  cleaved_par4 ,
            data = dat, subset = Container == "CBS Straw White")
summary(res)$coefficients

# green
res = coxph(Surv(days_os, PatientStatus_Description) ~ cleaved_par4, 
            data = dat, subset = Container == "CBS Straw Green")
summary(res)$coefficients












