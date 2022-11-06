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

# step 1. Get the rows with "doce" or "taxotere" (case insensitive) in the "txt_chemotherapy" column of DBBR.pheno.xlsx. 
# get the dates of the corresponding drug from the column.
filepath.pheno <-"/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/inputs/select_patients/DBBR.pheno.xlsx"
pheno.dat <- read_excel(path = filepath.pheno)

keywords <- c("doce" ,"taxotere")
txt_chemotherapy <- pheno.dat$Txt_Chemotherapy

match.info <- sapply(keywords, function(x) grepl(paste0(x), txt_chemotherapy, ignore.case = TRUE))
#match.info <- sapply(keywords, function(x) grepl(paste0("\\b", x, "\\b"), txt_chemotherapy, ignore.case = TRUE))
match_flag <- rowSums(match.info)
table(match_flag)
# match_flag
# 0   1 
# 404  63 

pheno.extract <- pheno.dat[match_flag> 0 , c("SubjectID", "Txt_Chemotherapy")]
dim(pheno.extract)


pheno.extract$eff.date <- rep(NA, dim(pheno.extract)[1])
pheno.extract$semi_colon <- rep(NA, dim(pheno.extract)[1])

for (row in 1:dim(pheno.extract)[1]) {
  strings <- unlist(strsplit(pheno.extract$Txt_Chemotherapy[row],split = " "))
  # to find if there are any dates in the current strings
  flag <- grepl(paste0("\\d{1,2}/\\d{1,2}/\\d{2,4}"), strings)
  
  word_1 <- grepl(keywords[1], strings, ignore.case = TRUE)
  word_2 <- grepl(keywords[2], strings, ignore.case = TRUE)
  chr_semicolon <- grepl(";", strings, ignore.case = TRUE)
  
  # if multiple keywords
  if ((sum(word_1) + sum(word_2) > 1) | (sum(flag) == 0) ) {
    next;
  }

  # no semicolon
  if (sum(chr_semicolon) == 0) {
    
    if (sum(flag) == 1) {
      pheno.extract$eff.date[row] <- str_extract(strings[flag], "[0-9]{1,2}/[0-9]{1,2}/[0-9]{2,4}")
    } else if (sum(flag) > 1) {
      # find positions of dates in current strings
      pos_dates <- which(flag)
      
      if (sum(word_1)>0) {
        this.word <- word_1
      } else if (sum(word_2)>0) {
        this.word <- word_2
      } 
      
      pos_word <- which(this.word)
      min_dis <- Inf; out.date <- NA
      for (pos_date in pos_dates) {
        if (abs(pos_date - pos_word) < min_dis  ) {
          min_dis <- abs(pos_date - pos_word)
          out.date <- strings[pos_date]
        } 
      }
      pheno.extract$eff.date[row] <- str_extract(out.date, "[0-9]{1,2}/[0-9]{1,2}/[0-9]{2,4}")
    }
  }
  
  # has semicolon
  if (sum(chr_semicolon) > 0) {
    
    pheno.extract$semi_colon[row] <- "YES"
    # then split by semicolon
    string_by_semicolon <- unlist(strsplit(pheno.extract$Txt_Chemotherapy[row],split = ";"))
    
    for (sub_string in string_by_semicolon) {
      
      sub.strings <- unlist(strsplit(sub_string, split = " "))

      sub.word_1 <- grepl(keywords[1], sub.strings, ignore.case = TRUE)
      sub.word_2 <- grepl(keywords[2], sub.strings, ignore.case = TRUE)
     
      if (sum(sub.word_1)>0) {
        sub.word <- sub.word_1
      } else if (sum(sub.word_2)>0) {
        sub.word <- sub.word_2
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
      pheno.extract$eff.date[row] <- str_extract(out.date, "[0-9]{1,2}/[0-9]{1,2}/[0-9]{2,4}")
    }
    
  }
  
}

WriteXLS::WriteXLS(pheno.extract,  ExcelFileName = "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/output_docetaxel/phenoFile.txt_chemotherapy.docetaxel.xls") 

#---------------------------------------------------------------------------------------------------------------------------------#
# further processed
filepath.pheno.out <-"/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/output_docetaxel/phenoFile.txt_chemotherapy.docetaxel.xls"
pheno.extract.dat <- read_excel(path = filepath.pheno.out)
pheno.extract.processed.dat <- pheno.extract.dat[- which(pheno.extract.dat$eff.date == "no date"), ]


pheno.extract.processed.dat$old_date <- pheno.extract.processed.dat$eff.date
#date.dat$drug_date_new <- rep(NA, dim(date.dat)[1])
num.date <- as.numeric(pheno.extract.processed.dat$eff.date)
num.date.index <- which(!is.na(num.date))
num_2_date <- as.Date(num.date, origin = "1899-12-30")
pheno.extract.processed.dat$eff.date[num.date.index] <- format(num_2_date[num.date.index], "%Y-%m-%d")

library(lubridate)
mdy <- mdy(pheno.extract.processed.dat$eff.date) 
ymd <- ymd(pheno.extract.processed.dat$eff.date)
mdy[is.na(mdy)] <- ymd[is.na(mdy)] # some dates are ambiguous, here we give 
pheno.extract.processed.dat$eff.date <- mdy        # mdy precedence over dmy

#sum(is.na(pheno.extract.processed.dat$eff.date))
# ID + date
pheno.extract.processed.dat <- pheno.extract.processed.dat %>% select(SubjectID, eff.date) 
WriteXLS::WriteXLS(pheno.extract.processed.dat,  ExcelFileName = "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/output_docetaxel/phenoFile.docetaxel.ID_Date.xls") 


##########################################################################################################################################
# step 2. Get the rows with "doce" or "taxotere" (case insensitive) in the "Treatment Text" column of qrySongTreatSeq.xlsx. 
# get the dates of the corresponding drug from the column. the entries are separated by ";" 
# and the date in each entry is usually in the front. for the rows do not have drug date in the "Treatment Text" column, 
# use the "DtRxStart" column

filepath.seq <-"/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/inputs/select_patients/qrySongTreatSeq.xlsx"
seq.dat <- read_excel(path = filepath.seq)

keywords <- c("doce" ,"taxotere")
trt <- seq.dat$`Treatment Text`

match.info <- sapply(keywords, function(x) grepl(paste0(x), trt, ignore.case = TRUE))
#match.info <- sapply(keywords, function(x) grepl(paste0("\\b", x, "\\b"), trt, ignore.case = TRUE))
match_flag <- rowSums(match.info)
table(match_flag)
# match_flag
# 0    1 
# 3142   72 

extract.dat <- seq.dat[match_flag> 0 , c("SUBJECTID", "Treatment Text", "DtRxStart")]
dim(extract.dat)

extract.dat$eff.date <- rep(NA, dim(extract.dat)[1])
extract.dat$semi_colon <- rep(NA, dim(extract.dat)[1])
for (row in 1:dim(extract.dat)[1]) {
  strings <- unlist(strsplit(extract.dat$`Treatment Text`[row],split = " "))
  # to find if there are any dates in the current strings
  flag <- grepl(paste0("\\d{1,2}/\\d{1,2}/\\d{1,2}"), strings)
  
  word_1 <- grepl(keywords[1], strings, ignore.case = TRUE)
  word_2 <- grepl(keywords[2], strings, ignore.case = TRUE)
  chr_semicolon <- grepl(";", strings, ignore.case = TRUE)
  
  # if multiple keywords
  if (sum(word_1) + sum(word_2) > 1) {
    next;
    # we will mannual fill the date
  }
  
  # if no date
  if (sum(flag) == 0) {
    # use the column "DtRxStart"
    extract.dat$eff.date[row] <- extract.dat$DtRxStart[row]
    
  }
  
  # no semicolon
  if (sum(chr_semicolon) == 0) {
    
    if (sum(flag) == 1) {
      extract.dat$eff.date[row] <- str_extract(strings[flag], "[0-9]{1,2}/[0-9]{1,2}/[0-9]{2,4}")
    } else if (sum(flag) > 1) {
      # find positions of dates in current strings
      pos_dates <- which(flag)
      
      if (sum(word_1)>0) {
        this.word <- word_1
      } else if (sum(word_2)>0) {
        this.word <- word_2
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
    string_by_semicolon <- unlist(strsplit(extract.dat$`Treatment Text`[row],split = ";"))
    
    for (sub_string in string_by_semicolon) {
      
      sub.strings <- unlist(strsplit(sub_string, split = " "))
      
      sub.word_1 <- grepl(keywords[1], sub.strings, ignore.case = TRUE)
      sub.word_2 <- grepl(keywords[2], sub.strings, ignore.case = TRUE)
      
      if (sum(sub.word_1)>0) {
        sub.word <- sub.word_1
      } else if (sum(sub.word_2)>0) {
        sub.word <- sub.word_2
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

WriteXLS::WriteXLS(extract.dat,  ExcelFileName = "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/output_docetaxel/seq.Treatment_Text.docetaxel.xls") 

#--------------------------------------------------------------------------------------------------------------------------------#
filepath.seq.out <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/output_docetaxel/seq.Treatment_Text.docetaxel.xls"
seq.extract.dat <- read_excel(path = filepath.seq.out)
seq.extract.processed.dat <- seq.extract.dat[-which(seq.extract.dat$eff.date == "no date") , ]


seq.extract.processed.dat$old_date <- seq.extract.processed.dat$eff.date
#date.dat$drug_date_new <- rep(NA, dim(date.dat)[1])
num.date <- as.numeric(seq.extract.processed.dat$eff.date)
num.date.index <- which(!is.na(num.date))
num_2_date <- as.Date(num.date, origin = "1899-12-30")
seq.extract.processed.dat$eff.date[num.date.index] <- format(num_2_date[num.date.index], "%Y-%m-%d")

library(lubridate)
mdy <- mdy(seq.extract.processed.dat$eff.date) 
ymd <- ymd(seq.extract.processed.dat$eff.date)
mdy[is.na(mdy)] <- ymd[is.na(mdy)] # some dates are ambiguous, here we give 
seq.extract.processed.dat$eff.date <- mdy        # mdy precedence over dmy

seq.extract.processed.dat <- seq.extract.processed.dat %>% select(SUBJECTID, eff.date) 
WriteXLS::WriteXLS(seq.extract.processed.dat,  ExcelFileName = "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/output_docetaxel/seq.docetaxel.ID_Date.xls") 


#######################################################################################################################################################
# step 3. put the dates and "SubjectID" you got in step1 and step2 together and sort the dates by "SubjectID", 
# use the earliest date for each patient

filepath.seq.ID_Date <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/output_docetaxel/seq.docetaxel.ID_Date.xls"
filepath.pheno.ID_Date <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/output_docetaxel/phenoFile.docetaxel.ID_Date.xls"
seq.ID_date <- read_excel(path = filepath.seq.ID_Date)
pheno.ID_date <- read_excel(path = filepath.pheno.ID_Date)
colnames(seq.ID_date)[1] <- "SubjectID"
stacked.ID_date <- rbind(seq.ID_date, pheno.ID_date)

stacked.earlist.dat <- stacked.ID_date %>% 
  group_by(SubjectID) %>% 
  arrange(eff.date) %>% slice(1L) 

#######################################################################################################################################################
# step 4. select "DOCETaxel" in "Order Name" column from "EHR___Orders" sheet of DEID_BDR 124420_Treatment_20210326.xlsx. 
# sort by "SubjectID",  then for the same "SubjectID" sort by the dates. take the earliest date for each "SubjectID"

filepath.sup <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/inputs/select_patients/DEID_BDR 124420_Treatment_20210326.xlsx"
sup.dat <- read_excel(path = filepath.sup, sheet = "EHR___Orders")
key.w <- "DOCETaxel"
trt <- sup.dat$`Order Name`

match.info <- sapply(key.w, function(x) grepl(paste0(x), trt, ignore.case = TRUE))
#match.info <- sapply(keywords, function(x) grepl(paste0("\\b", x, "\\b"), trt, ignore.case = TRUE))
match_flag <- rowSums(match.info)
table(match_flag)
# match_flag
# 0     1 
# 12763   278 

sup.dat <- sup.dat[match_flag>0, ]
sup.earlist.dat <- sup.dat %>% 
  group_by(SubjectID) %>% 
  arrange(`Date Task`) %>% slice(1L) %>% select(SubjectID, `Date Task`)

#######################################################################################################################################################
# step5. merge the dates you get from step3 with the dates you get from step4 and compare the consistency


full_joined.dat <- full_join(x=stacked.earlist.dat, y=sup.earlist.dat, by="SubjectID")
colnames(full_joined.dat)[2:3] <- c("Date_step3", "Date_step4")

full_joined.dat$joint_info <- NA
full_joined.dat$warning_info <- NA
for (i in 1:dim(full_joined.dat)[1]) {
  if (!is.na(full_joined.dat$Date_step3[i]) & !is.na(full_joined.dat$Date_step4[i])) {
    full_joined.dat$joint_info[i] <- "joint"
    
    if (full_joined.dat$Date_step3[i] == as.Date( full_joined.dat$Date_step4[i],'%m/%d/%Y')) {
      full_joined.dat$warning_info[i] <- "Equal"
    } else {
      full_joined.dat$warning_info[i] <- "Not Equal"
    }
  } else if (!is.na(full_joined.dat$Date_step3[i]) & is.na(full_joined.dat$Date_step4[i])) {
    full_joined.dat$joint_info[i] <- "only_pheno_seq"
  } else if (is.na(full_joined.dat$Date_step3[i]) & !is.na(full_joined.dat$Date_step4[i])) {
    full_joined.dat$joint_info[i] <- "only_sup"
  } else {
    full_joined.dat$joint_info[i] <- "both NA"
  }
}
WriteXLS::WriteXLS(full_joined.dat,  ExcelFileName = "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/output_docetaxel/full_join.docetaxel.xls") 


#-----------------------------------------------------------------------------------------------------#


filepath.joined <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/output_docetaxel/full_join.docetaxel.xls"
full_joined.dat <- read_excel(path = filepath.joined)

full_joined.dat$out_date <- NA
for (i in 1:dim(full_joined.dat)[1]) {
  if (full_joined.dat$joint_info[i] == "joint") {
    full_joined.dat$out_date[i] <- min(full_joined.dat$Date_pheno_seq[i], full_joined.dat$Date_sup[i])
  }
  if (full_joined.dat$joint_info[i] == "only_pheno_seq") {
    full_joined.dat$out_date[i] <- full_joined.dat$Date_pheno_seq[i]
  }
  if (full_joined.dat$joint_info[i] == "only_sup") {
    full_joined.dat$out_date[i] <- full_joined.dat$Date_sup[i]
  }
}

full_joined.dat <- full_joined.dat %>% select(SubjectID, out_date)

WriteXLS::WriteXLS(full_joined.dat,  ExcelFileName = "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/output_docetaxel/joined.docetaxel.ID_Date.xls") 

write.table(full_joined.dat, file = "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/output_docetaxel/joined.docetaxel.ID_Date.txt",
            row.names = FALSE, quote = FALSE) 
