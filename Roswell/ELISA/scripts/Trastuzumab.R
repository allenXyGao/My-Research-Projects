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


# DEID_BDR 124420_Treatment_20210326  sheet = "EHR___Orders"
filepath <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/inputs/select_patients/DEID_BDR 124420_Treatment_20210326.xlsx"
dat <- read_excel(path = filepath, sheet = "EHR___Orders")
keywords <- c("Ado-Trastuzumab", "Fam-Trastuzumab Deruxtecan", "Trastuzumab", "Trastuzumab (NCG195511)", "Pertuzumab")
dat.drug <- dat[dat$`Order Name` %in% keywords , ]

sum(is.na(dat.drug$`Date Task`))
dat.drug_ID_Date <- dat.drug %>% select(SubjectID, `Date Task`)


# DEID_BDR 124420_Treatment_20210326  Rx_Lapatinib
dat.Rx_Lapatinib <- read_excel(path = filepath, sheet = "Rx_Lapatinib")
dat.Rx_Lapatinib$tmp.date <-  as.Date(dat.Rx_Lapatinib$StartDate,'%m/%d/%Y')
dat.Rx_Lapatinib$year.StartDate <- as.numeric(format(dat.Rx_Lapatinib$tmp.date,'%Y'))

dat.Rx_Lapatinib <- dat.Rx_Lapatinib %>% filter(year.StartDate != 1900) %>% select(SubjectID, StartDate)


# add
colnames(dat.Rx_Lapatinib) <- c("SubjectID", "Date")
colnames(dat.drug_ID_Date) <- c("SubjectID", "Date")
dat.stacked <- rbind(dat.Rx_Lapatinib, dat.drug_ID_Date)

dat.stacked.earlist <- dat.stacked %>% 
  group_by(SubjectID) %>% 
  arrange(Date) %>% slice(1L)

# DBBR.pheno
filepath.pheno <-"/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/inputs/select_patients/DBBR.pheno.xlsx"
dat.pheno <- read_excel(path = filepath.pheno)

num.date <- as.numeric(dat.pheno$Anti.Her2.Date)
num.date.index <- which(!is.na(num.date))
num_2_date <- as.Date(num.date, origin = "1899-12-30")
dat.pheno$format.Anti.Her2.Date <- NA
dat.pheno$format.Anti.Her2.Date[num.date.index] <- format(num_2_date[num.date.index], "%Y-%m-%d")
dat.pheno.ID_DATE <- dat.pheno %>% select(SubjectID, format.Anti.Her2.Date)


# full join
full_joined.dat <- full_join(x=dat.pheno.ID_DATE, y=dat.stacked.earlist, by="SubjectID")

full_joined.dat$joint_info <- NA
full_joined.dat$warning_info <- NA
for (i in 1:dim(full_joined.dat)[1]) {
  if (!is.na(full_joined.dat$format.Anti.Her2.Date[i]) & !is.na(full_joined.dat$Date[i])) {
    full_joined.dat$joint_info[i] <- "joint"
    
    if (full_joined.dat$format.Anti.Her2.Date[i] == as.Date( full_joined.dat$Date[i],'%m/%d/%Y')) {
      full_joined.dat$warning_info[i] <- "Equal"
    } else {
      full_joined.dat$warning_info[i] <- "Not Equal"
    }
  } else if (!is.na(full_joined.dat$format.Anti.Her2.Date[i]) & is.na(full_joined.dat$Date[i])) {
    full_joined.dat$joint_info[i] <- "only pheno"
  } else if (is.na(full_joined.dat$format.Anti.Her2.Date[i]) & !is.na(full_joined.dat$Date[i])) {
    full_joined.dat$joint_info[i] <- "only sup"
  } else {
    full_joined.dat$joint_info[i] <- "both NA"
  }
}

WriteXLS::WriteXLS(full_joined.dat,  ExcelFileName = "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/output_Trastuzumab/Trastuzumab.full_join.xls") 
  

#----------------------------------------------------------------------------------------------------#
# output
filepath <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/output_Trastuzumab"
joined.dat <- read_excel(path = paste0(filepath, "/Trastuzumab.full_join.xls"))
joined.dat$Trastuzumab.Date <- NA

for (i in 1:dim(joined.dat)[1]) {
  if (!is.na(joined.dat$Anti.Her2.Date[i]) & !is.na(joined.dat$Date[i])) {
    joined.dat$Trastuzumab.Date[i] <- min(joined.dat$Anti.Her2.Date[i], joined.dat$Date[i])
  } else if (!is.na(joined.dat$Anti.Her2.Date[i]) & is.na(joined.dat$Date[i])) {
    joined.dat$Trastuzumab.Date[i] <- joined.dat$Anti.Her2.Date[i]
  } else if (is.na(joined.dat$Anti.Her2.Date[i]) & !is.na(joined.dat$Date[i])) {
    joined.dat$Trastuzumab.Date[i] <- joined.dat$Date[i]
  }
}


out_ID.Date <- joined.dat %>% select(SubjectID,Trastuzumab.Date) %>% filter(!is.na(Trastuzumab.Date))
WriteXLS::WriteXLS(out_ID.Date,  ExcelFileName = "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/output_Trastuzumab/Trastuzumab.ID_Date.xls") 




























