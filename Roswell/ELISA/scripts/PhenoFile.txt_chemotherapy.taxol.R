
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


filepath.pheno <-"/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/inputs/select_patients/DBBR.pheno.xlsx"
pheno.dat <- read_excel(path = filepath.pheno)

keywords <- c("paclitaxel", "taxol", "abraxane")
txt_chemotherapy <- pheno.dat$Txt_Chemotherapy

match.info <- sapply(keywords, function(x) grepl(paste0("\\b", x, "\\b"), txt_chemotherapy, ignore.case = TRUE))
match_flag <- rowSums(match.info)
table(match_flag)
# 0   1   2 
# 119 330  18 

extract.dat <- pheno.dat[match_flag> 0 , c("SubjectID", "Txt_Chemotherapy")]
dim(extract.dat)

extract.dat$eff.date <- rep(NA, dim(extract.dat)[1])
extract.dat$semi_colon <- rep(NA, dim(extract.dat)[1])
not.fill.rows <- c()
for (row in 1:dim(extract.dat)[1]) {
  strings <- unlist(strsplit(extract.dat$Txt_Chemotherapy[row],split = " "))
  # to find if there are any dates in the current strings
  flag <- grepl(paste0("\\d{1,2}/\\d{1,2}/\\d{1,2}"), strings)
  
  word_paclitaxel <- grepl("paclitaxel", strings, ignore.case = TRUE)
  word_taxol <- grepl("taxol", strings, ignore.case = TRUE)
  word_abraxane <- grepl("Abraxane", strings, ignore.case = TRUE)
  chr_semicolon <- grepl(";", strings, ignore.case = TRUE)
  
  # if multiple keywords
  if (sum(word_paclitaxel) + sum(word_taxol) + sum(word_abraxane) > 1) {
    not.fill.rows <- c(not.fill.rows, row)
    next;
    # we will mannual fill the date
  }
  
  # if no date
  if (sum(flag) == 0) {
    # use the column "DtRxStart"
    next;
    
  }
  
  # no semicolon
  if (sum(chr_semicolon) == 0) {
    
    if (sum(flag) == 1) {
      extract.dat$eff.date[row] <- str_extract(strings[flag], "[0-9]{1,2}/[0-9]{1,2}/[0-9]{2,4}")
    } else if (sum(flag) > 1) {
      # find positions of dates in current strings
      pos_dates <- which(flag)
      
      if (sum(word_abraxane)>0) {
        this.word <- word_abraxane
      } else if (sum(word_paclitaxel)>0) {
        this.word <- word_paclitaxel
      } else if (sum(word_taxol) > 0) {
        this.word <- word_taxol
      }
      
      pos_word <- which(this.word)
      min_dis <- Inf; out.date <- NA
      for (pos_date in pos_dates) {
        if (abs(pos_date - pos_word) < min_dis  ) {
          min_dis <- abs(pos_date - pos_word)
          out.date <- strings[pos_date]
        } 
      }
      extract.dat$eff.date[row] <- str_extract(out.date, "[0-9]{1,2}/[0-9]{1,2}/[0-9]{2,4}")
    }
  }

  # has semicolon
  if (sum(chr_semicolon) > 0) {
    
    extract.dat$semi_colon[row] <- "YES"
    # then split by semicolon
    string_by_semicolon <- unlist(strsplit(extract.dat$Txt_Chemotherapy[row],split = ";"))
    
    for (sub_string in string_by_semicolon) {
      
      sub.strings <- unlist(strsplit(sub_string, split = " "))
      
      sub.word_paclitaxel <- grepl("paclitaxel", sub.strings, ignore.case = TRUE)
      sub.word_taxol <- grepl("taxol", sub.strings, ignore.case = TRUE)
      sub.word_abraxane <- grepl("Abraxane", sub.strings, ignore.case = TRUE)
      
      if (sum(sub.word_abraxane)>0) {
        sub.word <- sub.word_abraxane
      } else if (sum(sub.word_paclitaxel)>0) {
        sub.word <- sub.word_paclitaxel
      } else if (sum(sub.word_taxol) > 0) {
        sub.word <- sub.word_taxol
      } else {
        next;
      }
      
      pos_dates <- which(grepl(paste0("\\d{1,2}/\\d{1,2}/\\d{1,2}"), sub.strings))
      pos_word <- which(sub.word)
      min_dis <- Inf; out.date <- NA
      for (pos_date in pos_dates) {
        if (abs(pos_date - pos_word) < min_dis  ) {
          min_dis <- abs(pos_date - pos_word)
          out.date <- sub.strings[pos_date]
        }
      }
      extract.dat$eff.date[row] <- str_extract(out.date, "[0-9]{1,2}/[0-9]{1,2}/[0-9]{2,4}")
    }
    
  }
    

  
}




WriteXLS::WriteXLS(extract.dat,  ExcelFileName = "/user/xinyugao/ELISA_project/phenoFile.txt_chemotherapy.taxol.xls") 



length(unique(extract.dat$SubjectID))




#######################################################################################

filepath.out <- "/user/xinyugao/ELISA_project/phenoFile.txt_chemotherapy.taxol.xls"
pheno.out <- read_excel(path = filepath.out)


pheno.out$old_date <- pheno.out$eff.date
#date.dat$drug_date_new <- rep(NA, dim(date.dat)[1])
num.date <- as.numeric(pheno.out$eff.date)
num.date.index <- which(!is.na(num.date))
num_2_date <- as.Date(num.date, origin = "1899-12-30")
pheno.out$eff.date[num.date.index] <- format(num_2_date[num.date.index], "%Y-%m-%d")

library(lubridate)
mdy <- mdy(pheno.out$eff.date) 
ymd <- ymd(pheno.out$eff.date)
mdy[is.na(mdy)] <- ymd[is.na(mdy)] # some dates are ambiguous, here we give 
pheno.out$eff.date <- mdy        # mdy precedence over dmy

sum(is.na(pheno.out$eff.date))
pheno.out <- pheno.out %>% select(-semi_colon)
colnames(pheno.out)[3:4] <- c("pheno.out_date", "pheno.old_date")

#########################################################################

filepath.taxol.out <- "/projects/rpci/qzhu/taxol_output/out.ID_Date.xls"
taxol.out <- read_excel(path = filepath.taxol.out)

colnames(taxol.out)
colnames(pheno.out)

full_joined.taxol_pheno.dat <- full_join(x=taxol.out, y=pheno.out, by=c("SUBJECTID"="SubjectID"))
full_joined.taxol_pheno.dat$joint_info <- NA
full_joined.taxol_pheno.dat$warning_info <- NA
for (i in 1:dim(full_joined.taxol_pheno.dat)[1]) {
  if (!is.na(full_joined.taxol_pheno.dat$out_date[i]) & !is.na(full_joined.taxol_pheno.dat$pheno.out_date[i])) {
    full_joined.taxol_pheno.dat$joint_info[i] <- "joint"
    if (full_joined.taxol_pheno.dat$out_date[i] != full_joined.taxol_pheno.dat$pheno.out_date[i]) {
      full_joined.taxol_pheno.dat$warning_info[i] <- "Not Equal"
    }
    if (full_joined.taxol_pheno.dat$out_date[i] == full_joined.taxol_pheno.dat$pheno.out_date[i]) {
      full_joined.taxol_pheno.dat$warning_info[i] <- "Equal"
    }
    
  } else if (!is.na(full_joined.taxol_pheno.dat$out_date[i]) & is.na(full_joined.taxol_pheno.dat$pheno.out_date[i])) {
    full_joined.taxol_pheno.dat$joint_info[i] <- "taxol_output_only"
  } else if (is.na(full_joined.taxol_pheno.dat$out_date[i]) & !is.na(full_joined.taxol_pheno.dat$pheno.out_date[i])) {
    full_joined.taxol_pheno.dat$joint_info[i] <- "pheno_out_only"
  }
}



WriteXLS::WriteXLS(full_joined.taxol_pheno.dat,  ExcelFileName = "/user/xinyugao/ELISA_project/full_joined.taxol_pheno.dat.xls") 

################################################################################################################################################

filepath.joined.pheno_taxol.out <-"/user/xinyugao/ELISA_project/full_joined.taxol_pheno.dat.xls"
dat.joined.pheno_taxol.out  <- read_excel(path = filepath.joined.pheno_taxol.out)
dat.joined.pheno_taxol.out <- dat.joined.pheno_taxol.out[-388, ] 

dat.joined.pheno_taxol.out$taxol_updated_out_date <- dat.joined.pheno_taxol.out$taxol_out_date
#dat.neq <- dat.joined.pheno_taxol.out[!is.na(dat.joined.pheno_taxol.out$warning_info) & dat.joined.pheno_taxol.out$warning_info == "Not Equal", ]
dat.joined.pheno_taxol.out$taxol_updated_out_date[dat.joined.pheno_taxol.out$SUBJECTID == 'PT-00232610'] <- 
  dat.joined.pheno_taxol.out$pheno.out_date[dat.joined.pheno_taxol.out$SUBJECTID == 'PT-00232610']

dat.joined.pheno_taxol.out$taxol_updated_out_date[dat.joined.pheno_taxol.out$SUBJECTID == 'PT-00310352'] <- 
  dat.joined.pheno_taxol.out$pheno.out_date[dat.joined.pheno_taxol.out$SUBJECTID == 'PT-00310352']

dat.joined.pheno_taxol.out$taxol_updated_out_date[dat.joined.pheno_taxol.out$SUBJECTID == 'PT-00344607'] <-
  dat.joined.pheno_taxol.out$pheno.out_date[dat.joined.pheno_taxol.out$SUBJECTID == 'PT-00344607']


dat.joined.pheno_taxol.out$taxol_updated_out_date[dat.joined.pheno_taxol.out$SUBJECTID == 'PT-00183817'] <- 
  dat.joined.pheno_taxol.out$pheno.out_date[dat.joined.pheno_taxol.out$SUBJECTID == 'PT-00183817']
dat.joined.pheno_taxol.out$taxol_updated_out_date[dat.joined.pheno_taxol.out$SUBJECTID == 'PT-00137254'] <- 
  dat.joined.pheno_taxol.out$pheno.out_date[dat.joined.pheno_taxol.out$SUBJECTID == 'PT-00137254']


WriteXLS::WriteXLS(dat.joined.pheno_taxol.out,  ExcelFileName = "/user/xinyugao/ELISA_project/full_joined.taxol_pheno.dat.xls") 


taxol.ID_Date_only.updated_Out <- dat.joined.pheno_taxol.out %>% select(SUBJECTID, taxol_updated_out_date)
WriteXLS::WriteXLS(taxol.ID_Date_only.updated_Out ,  ExcelFileName = "/user/xinyugao/ELISA_project/updated.Taxol.ID_Date.only.xls") 
