# 1. select "DOXOrubicin", "DOXOrubicin Liposome", "DOXOrubicin Liposome.", "Epirubicin" in "Order Name" column.
# 2. sort by "SubjectID", then for the same "SubjectID" sort by "Date Task"
# 3. take the earliest "Date Task" for each "SubjectID"
# 4. merge the results with DBBR.pheno.xlsx and compare the consistency with "Doxo.Date" in this file.

rm(list=ls())
mypaths <- .libPaths()
mypaths <- c(mypaths, "/projects/rpci/qzhu/xinyu/Software/R/R-4.0.2")
.libPaths(mypaths)

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
#install.packages("tidytext", lib="/projects/rpci/qzhu/xinyu/Software/R/R-4.0.2")
library(tidytext, lib.loc = "/projects/rpci/qzhu/xinyu/Software/R/R-4.0.2")


# DEID_BDR 124420_Treatment_20210326
filepath <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/inputs/select_patients/DEID_BDR 124420_Treatment_20210326.xlsx"
dat <- read_excel(path = filepath, sheet = "EHR___Orders")
keywords <- c( "DOXOrubicin", "DOXOrubicin Liposome", "DOXOrubicin Liposome.", "Epirubicin")
dat.drug <- dat[dat$`Order Name` %in% keywords , ]

sum(is.na(dat.drug$`Date Task`))
# 0

length(unique(dat.drug$SubjectID))
# 292 unique ID

earlist.dat <- dat.drug %>% 
  group_by(SubjectID) %>% 
  arrange(`Date Task`) %>% slice(1L) %>% select(SubjectID, `Date Task`)
#---------------------------------------------------------------------------------------------------------#


# DBBR.pheno
filepath.sup <-"/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/inputs/select_patients/DBBR.pheno.xlsx"
dat.sup <- read_excel(path = filepath.sup)
# dim 467 * 54
length(unique(dat.sup$SubjectID))
# 467 unique ID


# dat.sup$Txt_Chemotherapy[dat.sup$SubjectID == "PT-00165116"]
# # "3/31/09 Paclitaxel. Taxol, Herceptin, Xeloda, lapatinib, exemestane, Afinitor, doxorubicin, Eribulin and now Ixempra."
# dat.sup$Txt_Chemotherapy[dat.sup$SubjectID == "PT-00206983"]
# # "5 cycles of AC through 01/2010"
# dat.sup$Txt_Chemotherapy[dat.sup$SubjectID == "PT-00233927"]
# # "10/20/11(RPCI) TAXOL, AVASTIN,  4 cycles of Adriamycin and cytoxan on 3/2012"
# dat.sup$Txt_Chemotherapy[dat.sup$SubjectID == "PT-00185591"]
# # "1/3/2007 AC & TAXOL& HERCEPTIN (RPCI)"
# dat.sup$Txt_Chemotherapy[dat.sup$SubjectID == "PT-00091950"]
# # "3/19/10-RPCI-Adriamycin, Cytoxan,  then Taxol"
# dat.sup$Txt_Chemotherapy[dat.sup$SubjectID == "PT-00184210"]
# # "11/8/06 RPCI: ADRIAMYCIN, CYTOXAN; 1/16/07 RPCI: TAXOL"
# dat.sup$Txt_Chemotherapy[dat.sup$SubjectID == "PT-00186546"]
# # "3/12/07-RPCI-ADRIAMYCIN,CYTOXAN; 5/8/07-RPCI-TAXOL"
# dat.sup$Txt_Chemotherapy[dat.sup$SubjectID == "PT-00177895"]
# # "12/28/2005 AC & TAXOL (RPCI)"
# dat.sup$Txt_Chemotherapy[dat.sup$SubjectID == "PT-00344607"]
# # "1/99/17 (est.month)(Dr. Steinbrenner in VA) ddACT under care of , last cycle complicated by perianal abscess; 4/3/2017 dd Paclitaxel (neoadjuvant); 
# # 4/28/17 (RP) ddACT: Paclitaxel (neoadjuvant, transfer of care); 5/26/17 (RP) Abraxane (d/t side effects on"
# dat.sup$Txt_Chemotherapy[dat.sup$SubjectID == "PT-00260077"]
# # "06/05/2013 Adriamycin and Cytoxan"
# dat.sup$Txt_Chemotherapy[dat.sup$SubjectID == "PT-00114899"]
# # "9/15/11(RPCI) PACLITAXEL, TRASTUZUMAB"
# dat.sup$Txt_Chemotherapy[dat.sup$SubjectID == "PT-00375458"]

num.date <- as.numeric(dat.sup$Doxo.Date)
num.date.index <- which(!is.na(num.date))
num_2_date <- as.Date(num.date, origin = "1899-12-30")
dat.sup$Doxo.Date[num.date.index] <- format(num_2_date[num.date.index], "%Y-%m-%d")
dat.sup$Doxo.Date.missing <- "NO"
for (i in 1: dim(dat.sup)[1]) {
  if (dat.sup$Doxo.Date[i] == "NA") {
    dat.sup$Doxo.Date.missing[i] <- "YES"
  }
}
table(dat.sup$Doxo.Date.missing)
# NO YES 
# 383  84 


sup.earlist.dat <- dat.sup %>% 
  group_by(SubjectID) %>% 
  arrange(Doxo.Date) %>% slice(1L) %>% select(SubjectID, Doxo.Date, Doxo.Date.missing)
#---------------------------------------------------------------------------------------------------------#

## merge
full_joined.dat <- full_join(x=earlist.dat, y=sup.earlist.dat, by="SubjectID")
full_joined.dat$`Date Task` <- format(full_joined.dat$`Date Task`, "%Y-%m-%d")

full_joined.dat$joint_info <- NA
full_joined.dat$warning_info <- NA
for (i in 1:dim(full_joined.dat)[1]) {
  if (!is.na(full_joined.dat$`Date Task`[i]) & full_joined.dat$Doxo.Date.missing[i] == "NO") {
    full_joined.dat$joint_info[i] <- "joint"
    if (full_joined.dat$`Date Task`[i] == full_joined.dat$Doxo.Date[i]) {
      full_joined.dat$warning_info[i] <-"Equal"
    }
    if (full_joined.dat$`Date Task`[i] != full_joined.dat$Doxo.Date[i]) {
      full_joined.dat$warning_info[i] <-"Not Equal"
      #full_joined.dat$abs_diff_days[i] <- abs(full_joined.dat$`Date Task`[i] - full_joined.dat$Doxo.Date[i])
      }
  } else if (!is.na(full_joined.dat$`Date Task`[i]) & full_joined.dat$Doxo.Date.missing[i] == "NO") {
    full_joined.dat$joint_info[i] <- "joint"
    full_joined.dat$warning_info[i] <-"Joint, but Doxo.Date missing"
  } else if (!is.na(full_joined.dat$`Date Task`[i]) & full_joined.dat$Doxo.Date.missing[i] == "YES") {
    full_joined.dat$joint_info[i] <- "Only Date Task"
  } else if (is.na(full_joined.dat$`Date Task`[i]) & full_joined.dat$Doxo.Date.missing[i] == "NO") {
    full_joined.dat$joint_info[i] <- "Only Doxo.Date"
  }
}

full_joined.dat$Doxo.Date[full_joined.dat$SubjectID == "PT-00260077"] <- format(as.Date("2013-06-13", format="%Y-%m-%d"), "%Y-%m-%d")
#WriteXLS::WriteXLS(full_joined.dat,  ExcelFileName = "/user/xinyugao/ELISA_project/out.doxorubicin.full_join.xls") 




full_joined.dat$joined.date <- NA
for (i in 1:dim(full_joined.dat)[1]) {
  if (is.na(full_joined.dat$joint_info[i] )) {
    next;
  }
  if (full_joined.dat$joint_info[i] == "joint") {
    full_joined.dat$joined.date[i] <- min(full_joined.dat$`Date Task`[i], full_joined.dat$Doxo.Date[i])
  }
  if (full_joined.dat$joint_info[i] == "Only Date Task") {
    full_joined.dat$joined.date[i] <- full_joined.dat$`Date Task`[i]
  }
  if (full_joined.dat$joint_info[i] == "Only Doxo.Date") {
    full_joined.dat$joined.date[i] <- full_joined.dat$Doxo.Date[i]
  }
}

out <- full_joined.dat[-which(is.na(full_joined.dat$joined.date)),]

# filepath.taxol <-"/user/xinyugao/ELISA_project/SAMPLE/out.ID_Date.xls"
# dat.taxol <- read_excel(path = filepath.taxol)
# 
# rubicin.joint.neq <- out[out$joint_info == "joint" & out$warning_info == "Not Equal", ]
# 
# inner_join.dat <- inner_join(rubicin.joint.neq, dat.taxol, by=c("SubjectID"="SUBJECTID"))
# for (i in 1:dim(inner_join.dat)[1]) {
#   if (inner_join.dat$out_date[i] > min(inner_join.dat$`Date Task`[i], inner_join.dat$Doxo.Date[i]) &
#       inner_join.dat$out_date[i] < max(inner_join.dat$`Date Task`[i], inner_join.dat$Doxo.Date[i]))
#     print(inner_join.dat$SubjectID[i])
# }
colnames(out)[7] <- "doxorubicin.Date"
out.selected <- out %>% select("SubjectID", "doxorubicin.Date")
  
WriteXLS::WriteXLS(out,  ExcelFileName = "/user/xinyugao/ELISA_project/out.doxorubicin.full_join.processed.xls") 
WriteXLS::WriteXLS(out.selected,  ExcelFileName = "/user/xinyugao/ELISA_project/out.doxorubicin.full_join.ID_Date.xls") 




