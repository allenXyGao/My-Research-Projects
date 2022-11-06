# 1. Find the rows with "paclitaxel" or "taxol" (case insensitive) in the "Treatment Text" column. 
# 2. get the dates of the corresponding drug from the column. Please note that the date is usually right before the drug name.
#    However, there're some cases like "7/22/16 (RPCI) AC-T: Cyclophosphamide, Doxorubicin, followed by Taxol 9/16/16" and 
#    "7/9/07 - RPCI - DOXORUBICIN & CYCLOPHOSPHAMIDE, FOLLOWED BY TAXOL (9/3/07)", where the date is after the drug name.
# 3. for the rows do not have drug date in the "Treatment Text" column, use the "DtRxStart" column

# Your results should be a table with 2 columns: SUBJECTID and drug date. Some patients may have multiple drug date and therefore multiple rows.
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

filepath <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/inputs/select_patients/qrySongTreatSeq.xlsx"
treat.dat <- read_excel(path = filepath)
# step 1: find rows containing keywords (case insensitive)
keywords <- c("paclitaxel", "taxol", "Abraxane")
trt <- treat.dat$`Treatment Text`

match.info <- sapply(keywords, function(x) grepl(paste0("\\b", x, "\\b"), trt, ignore.case = TRUE))
match_flag <- rowSums(match.info)
table(match_flag)
match_flag
#   0    1    2    3 
# 2840  357   16    1 


extract.dat <- treat.dat[match_flag> 0 , c("SUBJECTID", "DtRxStart", "Treatment Text")]
dim(extract.dat)
# 374   3
# [1] "Taxol (Paclitaxel)"                                                                                                                                                                                   
# [2] "8/29/16 (RPCI) Paclitaxel (neoadjuvant) (D/t PT`s comorbidities and increased risk of pneumonia during the winter, aggressive chemo was not recommended. Taxol d/c`d 11/21/16 s/p 9/12 planned doses)"
# [3] "4/28/17 (RP) ddACT: Paclitaxel (neoadjuvant, transfer of care); 5/26/17 (RP) Abraxane (d/t side effects on Taxol); 11/3/17 (RP) Carboplatin (adjuvant, did not get cycle 4 d/t thrombocytopenia)"

# step 2:get the dates of the corresponding drug from the column
extract.dat$eff.date <- rep(NA, dim(extract.dat)[1])
extract.dat$semi_colon <- rep(NA, dim(extract.dat)[1])
not.fill.rows <- c()
for (row in 1:dim(extract.dat)[1]) {
  strings <- unlist(strsplit(extract.dat$`Treatment Text`[row],split = " "))
  # to find if there are any dates in the current strings
  flag <- grepl(paste0("\\d{1,2}/\\d{1,2}/\\d{1,2}"), strings)
  
  word_paclitaxel <- grepl("paclitaxel", strings, ignore.case = TRUE)
  word_taxol <- grepl("taxol", strings, ignore.case = TRUE)
  word_abraxane <- grepl("Abraxane", strings, ignore.case = TRUE)
  chr_semicolon <- grepl(";", strings, ignore.case = TRUE)
  
  if (sum(word_paclitaxel) + sum(word_taxol) + sum(word_abraxane) > 1) {
    not.fill.rows <- c(not.fill.rows, row)
    next;
    # we will mannual fill the date
  }
  
  
  # case1: no date
  if (sum(flag) == 0) {
    # use the column "DtRxStart"
    extract.dat$eff.date[row] <- extract.dat$DtRxStart[row]
    
  } # case 2: only 1 date
  else if (sum(flag) == 1) {
    
    if (sum(chr_semicolon) > 0) {
      extract.dat$semi_colon[row] <- "YES"
    }
    # use this date
    extract.dat$eff.date[row] <- str_extract(strings[flag], "[0-9]{1,2}/[0-9]{1,2}/[0-9]{1,2,3,4}")
  }  
  # case 3: multiple dates, select the correct date
  else if (sum(flag) > 1) {
    
    if (sum(chr_semicolon) > 0) {
      extract.dat$semi_colon[row] <- "YES"
    }
    
    # find positions of dates in current strings
    pos_dates <- which(flag)
    
    # paclitaxel
    if ( sum(word_paclitaxel) > 0 ) {
      pos_word <- which(word_paclitaxel)
      min_dis <- Inf; out.date <- NA
      # when abs(dis1) = abs(dis2), select the first one, becasue the for loop
      # starts from the begining of the string, so the first nearest date (before the word) is selected
      for (pos_date in pos_dates) {
        if (abs(pos_date - pos_word) < min_dis  ) {
          min_dis <- abs(pos_date - pos_word)
          out.date <- strings[pos_date]
        } 
      }
      
    } 
    # taxol
    if ( sum(word_taxol) > 0 ) {
      pos_word <- which(word_taxol)
      min_dis <- Inf; out.date <- NA
      for (pos_date in pos_dates) {
        if (abs(pos_date - pos_word) < min_dis  ) {
          min_dis <- abs(pos_date - pos_word)
          out.date <- strings[pos_date]
        }
      }
    }
    
    # abraxane
    if ( sum(word_abraxane) > 0 ) {
      pos_word <- which(word_abraxane)
      min_dis <- Inf; out.date <- NA
      for (pos_date in pos_dates) {
        if (abs(pos_date - pos_word) < min_dis  ) {
          min_dis <- abs(pos_date - pos_word)
          out.date <- strings[pos_date]
        }
      }
    }
    
    
    extract.dat$eff.date[row] <- str_extract(out.date, "[0-9]{1,2}/[0-9]{1,2}/[0-9]{1,2,3,4}")
  }
}


#str_extract(extract.dat$eff.date[5], "[0-9]{1,2}/[0-9]{1,2}/[0-9]{1,2}")

colnames(extract.dat)[4] <- "drug_date"
extract.dat$not.fill.rows <- rep(NA, dim(extract.dat)[1])
extract.dat$not.fill.rows[1:length(not.fill.rows)] <- not.fill.rows

#write.csv(extract.dat, file = "/user/xinyugao/ELISA_project/out.csv" )
WriteXLS::WriteXLS(extract.dat,  ExcelFileName = "/user/xinyugao/ELISA_project/out.3keys.xls", SheetNames = "main") 



# 
# 
# strings <- unlist(strsplit(extract.dat$`Treatment Text`[1],split = " "))
# flag <- grepl(paste0("\\d{1,2}/\\d{1,2}/\\d{1,2}"), strings)
# # case1: no date
# if (sum(flag) == 0) {
#   # use the column "DtRxStart"
# } else if (sum(flag) == 1) {
#   # use this date
# }  else if (sum(flag) > 1) {
#   # find positions of dates in current strings
#   pos_dates <- which(flag)
#   
#   word_paclitaxel <- grepl("paclitaxel", strings, ignore.case = TRUE)
#   word_taxol <- grepl("taxol", strings, ignore.case = TRUE)
#   # here we need to pre-exclude paclitaxel and taxol appear simutaneously
#   # paclitaxel
#   if ( sum(word_paclitaxel) > 0 ) {
#     pos_word <- which(word_paclitaxel)
#     min_dis <- Inf; out.date <- NA
#     # when abs(dis1) = abs(dis2), select the first one, becasue the for loop
#     # starts from the begining of the string, so the first nearest date (before the word) is selected
#     for (pos_date in pos_dates) {
#       if (abs(pos_date - pos_word) < min_dis  ) {
#         min_dis <- abs(pos_date - pos_word)
#         out.date <- strings[pos_date]
#       } 
#     }
#     
#   } 
#   # taxol
#   if ( sum(word_taxol) > 0 ) {
#     pos_word <- which(word_taxol)
#     min_dis <- Inf; out.date <- NA
#     for (pos_date in pos_dates) {
#       if (abs(pos_date - pos_word) < min_dis  ) {
#         min_dis <- abs(pos_date - pos_word)
#         out.date <- strings[pos_date]
#       }
#     }
#   }
#   
#   
# }


#######################################################################################################################################################
# 1. select "PACLItaxel" and "PACLItaxel (Protein-Bound)" in "Order Name" column.
# 
# 2. sort by "SubjectID", then for the same "SubjectID" sort by "Date Task"
# 
# 3. take the earliest "Date Task" for each "SubjectID"
# 
# 4. merge the results with the results you got from qrySongTreatSeq.xlsx and compare the consistency of the dates. 
# 
# For the inconsistent incidences, in theory the dates from qrySongTreatSeq.xlsx should be earlier than from DEID_BDR 124420_Treatment_20210326.xlsx.

filepath.sup <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/inputs/select_patients/DEID_BDR 124420_Treatment_20210326.xlsx"
sup.dat <- read_excel(path = filepath.sup, sheet = "EHR___Orders")
sup.dat.drug <- sup.dat[sup.dat$`Order Name` == "PACLItaxel" | sup.dat$`Order Name` =="PACLItaxel (Protein-Bound)", ]

head(sup.dat.drug)

sup.earlist.dat <- sup.dat.drug %>% 
  group_by(SubjectID) %>% 
  arrange(`Date Task`) %>% slice(1L) %>% select(SubjectID, `Date Task`)


filepath.out <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/inputs/select_patients/out.3keys_new.xls"
date.dat <- read_excel(path = filepath.out, col_types = c("text", "text", "text", "text", "text"))


length(unique(date.dat$SUBJECTID))
# 348
# mark ACT and NO
ind_ACT <- which(date.dat$drug_date %in% c("ACT"))
ind_no <- which(date.dat$drug_date %in% c("no"))

date.dat$old_date <- date.dat$drug_date
#date.dat$drug_date_new <- rep(NA, dim(date.dat)[1])
num.date <- as.numeric(date.dat$drug_date)
num.date.index <- which(!is.na(num.date))
num_2_date <- as.Date(num.date, origin = "1899-12-30")
date.dat$drug_date[num.date.index] <- format(num_2_date[num.date.index], "%Y-%m-%d")

library(lubridate)
mdy <- mdy(date.dat$drug_date) 
ymd <- ymd(date.dat$drug_date)
mdy[is.na(mdy)] <- ymd[is.na(mdy)] # some dates are ambiguous, here we give 
date.dat$drug_date <- mdy        # mdy precedence over dmy

#date.dat$drug_date[ind_ACT] <- "2222-02-02"
date.dat$drug_date[ind_no] <- NA

ind.day99 <- which(is.na(date.dat$drug_date))
for (ind in ind.day99){
  cur_date <- unlist(strsplit(date.dat$old_date[ind],split = "/"))
  date.dat$drug_date[ind] <- as.Date(paste0(cur_date[1], "/30/20", cur_date[3]), "%m/%d/%Y")
}

sum(is.na(date.dat$drug_date))
      
# find the earlist record
date.earlist.dat <- date.dat %>% 
  group_by(SUBJECTID) %>% 
  arrange(drug_date) %>% slice(1L)


full_joined.dat <- full_join(x=date.earlist.dat, y=sup.earlist.dat, by=c("SUBJECTID"="SubjectID"))
full_joined.dat$joint_info <- NA
full_joined.dat$warning_info <- NA
for (i in 1:dim(full_joined.dat)[1]) {
  if (!is.na(full_joined.dat$drug_date[i]) & !is.na(full_joined.dat$`Date Task`[i])) {
    full_joined.dat$joint_info[i] <- "joint"
    if (full_joined.dat$drug_date[i] > full_joined.dat$`Date Task`[i]) {
      full_joined.dat$warning_info[i] <- "Warning"
    }
    if (full_joined.dat$drug_date[i] == full_joined.dat$`Date Task`[i]) {
      full_joined.dat$warning_info[i] <- "Equal"
    }
    
  } else if (!is.na(full_joined.dat$drug_date[i]) & is.na(full_joined.dat$`Date Task`[i])) {
    full_joined.dat$joint_info[i] <- "trt_only"
  } else if (is.na(full_joined.dat$drug_date[i]) & !is.na(full_joined.dat$`Date Task`[i])) {
    full_joined.dat$joint_info[i] <- "sup_only"
  }
}

full_joined.dat$drug_date[full_joined.dat$SUBJECTID =="PT-00193689"] <- "2008-02-28"
full_joined.dat$warning_info[full_joined.dat$SUBJECTID =="PT-00193689"] <- "Equal"

WriteXLS::WriteXLS(full_joined.dat,  ExcelFileName = "/user/xinyugao/ELISA_project/out.full_join.xls") 

####################################################################################################################


full_joined.dat <- full_joined.dat[- which(is.na(full_joined.dat$drug_date) & is.na(full_joined.dat$`Date Task`)), ]

full_joined.dat$out_date <- NA
for (i in 1:dim(full_joined.dat)[1]) {
  if (full_joined.dat$joint_info[i] == "joint") {
    full_joined.dat$out_date[i] <- format(min(full_joined.dat$drug_date[i], full_joined.dat$`Date Task`[i]), "%Y-%m-%d")
  }
  if (full_joined.dat$joint_info[i] == "trt_only") {
    full_joined.dat$out_date[i] <- format(full_joined.dat$drug_date[i], "%Y-%m-%d")
  }
  if (full_joined.dat$joint_info[i] == "sup_only") {
    full_joined.dat$out_date[i] <- format(full_joined.dat$`Date Task`[i], "%Y-%m-%d")
  }
}

out.dat <- full_joined.dat %>% select(SUBJECTID, out_date)

WriteXLS::WriteXLS(out.dat,  ExcelFileName = "/user/xinyugao/ELISA_project/out.ID_Date.xls") 
