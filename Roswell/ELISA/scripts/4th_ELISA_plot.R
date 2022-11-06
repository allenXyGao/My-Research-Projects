rm(list=ls())
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

# please draw plots for using 4thELISA_450nm-570nm, 4thELISA_450nm-540nm, 4thELISA_450nm-570nm.30K, 4thELISA_450nm-540nm.30K. 
# Sample information is in "Sheet1" sheet of 4th_ELISA_samples_results.xlsx.
# please connect the samples from the same patient in the plots like what you did for the plots of 1st run.


filepath.concentration.540nm <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/inputs/4th/4thELISA_450nm-540nm"
filepath.concentration.570nm <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/inputs/4th/4thELISA_450nm-570nm"
filepath.concentration.540nm.30K <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/inputs/4th/4thELISA_450nm-540nm.30K"
filepath.concentration.570nm.30K <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/inputs/4th/4thELISA_450nm-570nm.30K"

filepath.sample.info <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/inputs/4th/4th_ELISA_samples_results.xlsx"# sample info
#filepath.sample.info.1st <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/inputs/ELISA_samples_barcode.xlsx"

sample.info <- read_excel(path = filepath.sample.info, sheet = "Sheet1")
selected.columns <- c("SubjectID", "BARCODE", "TXHISTORYSTATUS" , "rs6494889")
sample.info.dat <- sample.info[, selected.columns]
dim(sample.info.dat)
# 16 * 4
#sample.info.1st <- read_excel(path = filepath.sample.info.1st, sheet = "qryPIyaoSamples (2)")
apply(sample.info.dat, 2, function(x) length(unique(x)))
# SubjectID         BARCODE TXHISTORYSTATUS       rs6494889 
# 9              16               3               3 
# table(sample.info.dat$BARCODE)[table(sample.info.dat$BARCODE) > 1]

table(sample.info.dat$TXHISTORYSTATUS)
# No prior treatment (excluding diagnostic biopsy) Post-surgical AND Post-adjuvant or systemic therapy       Post-surgical/No adjuvant or systemic therapy 
# 7                                                   5                                                   4 

sample.info.dat$TXHISTORYSTATUS[sample.info.dat$TXHISTORYSTATUS == "No prior treatment (excluding diagnostic biopsy)" | 
                                  sample.info.dat$TXHISTORYSTATUS == "NO PRIOR TREATMENT (EXCLUDING DIAGNOSTIC BIOPSY)"] <- "No Prior"
sample.info.dat$TXHISTORYSTATUS[sample.info.dat$TXHISTORYSTATUS == "Post-surgical/No adjuvant or systemic therapy"  |
                                  sample.info.dat$TXHISTORYSTATUS ==  "POST-SURGICAL/NO ADJUVANT OR SYSTEMIC THERAPY" |
                                  sample.info.dat$TXHISTORYSTATUS == "Post-surgical/No adjuvant or systemic therapy "] <- "Post but no-adjuvant"
sample.info.dat$TXHISTORYSTATUS[sample.info.dat$TXHISTORYSTATUS == "POST-SURGICAL AND POST-ADJUVANT OR SYSTEMIC THERAPY" | 
                                  sample.info.dat$TXHISTORYSTATUS == "Post-surgical AND Post-adjuvant or systemic therapy"] <- "Post & Post-adjuvant"
table(sample.info.dat$TXHISTORYSTATUS)
# No Prior Post but no-adjuvant Post & Post-adjuvant 
# 7                    4                    5 

sample.info.dat$variant.info <- rep(NA, dim(sample.info.dat)[1])
for (i in 1:dim(sample.info.dat)[1]) {
  if (is.na(sample.info.dat$rs6494889[i])) {sample.info.dat$variant.info[i] <- "freshly frozen"}
  if (!is.na(sample.info.dat$rs6494889[i]) & sample.info.dat$rs6494889[i] == 0) {sample.info.dat$variant.info[i] <- "no variant"}
  if (!is.na(sample.info.dat$rs6494889[i]) & (sample.info.dat$rs6494889[i] == 1 | sample.info.dat$rs6494889[i] == 2)) {sample.info.dat$variant.info[i] <- "has variant"}
}
table(sample.info.dat$variant.info)
# freshly frozen    has variant     no variant 
# 3                    1             12 


sample.info.dat$TXHISTORYSTATUS <- factor(sample.info.dat$TXHISTORYSTATUS, levels = c("No Prior", "Post but no-adjuvant", 
                                                                                      "Post & Post-adjuvant"))
sample.info.dat$variant.info <- factor(sample.info.dat$variant.info, levels = c("freshly frozen", "no variant", "has variant"))
sample.info.dat.easy_check <- sample.info.dat %>% arrange(TXHISTORYSTATUS, variant.info)

sample.info.dat.easy_check$x_numeric <- rep(NA, dim(sample.info.dat)[1])
for (i in 1:dim(sample.info.dat)[1]) {
  if (sample.info.dat.easy_check$TXHISTORYSTATUS[i] == "No Prior" & sample.info.dat.easy_check$variant.info[i] == "freshly frozen") {
    sample.info.dat.easy_check$x_numeric[i] <- 1
  }
  else if (sample.info.dat.easy_check$TXHISTORYSTATUS[i] == "No Prior" & sample.info.dat.easy_check$variant.info[i] == "no variant") {
    sample.info.dat.easy_check$x_numeric[i] <- 2
  }
  else if (sample.info.dat.easy_check$TXHISTORYSTATUS[i] == "No Prior" & sample.info.dat.easy_check$variant.info[i] == "has variant") {
    sample.info.dat.easy_check$x_numeric[i] <- 3
  }
  else if (sample.info.dat.easy_check$TXHISTORYSTATUS[i] == "Post but no-adjuvant" & sample.info.dat.easy_check$variant.info[i] == "freshly frozen") {
    sample.info.dat.easy_check$x_numeric[i] <- 4
  }
  else if (sample.info.dat.easy_check$TXHISTORYSTATUS[i] == "Post but no-adjuvant" & sample.info.dat.easy_check$variant.info[i] == "no variant") {
    sample.info.dat.easy_check$x_numeric[i] <- 5
  }
  else if (sample.info.dat.easy_check$TXHISTORYSTATUS[i] == "Post but no-adjuvant" & sample.info.dat.easy_check$variant.info[i] == "has variant") {
    sample.info.dat.easy_check$x_numeric[i] <- 6
  }
  else if (sample.info.dat.easy_check$TXHISTORYSTATUS[i] == "Post & Post-adjuvant" & sample.info.dat.easy_check$variant.info[i] == "freshly frozen") {
    sample.info.dat.easy_check$x_numeric[i] <- 7
  }
  else if (sample.info.dat.easy_check$TXHISTORYSTATUS[i] == "Post & Post-adjuvant" & sample.info.dat.easy_check$variant.info[i] == "no variant") {
    sample.info.dat.easy_check$x_numeric[i] <- 8
  }
  else if (sample.info.dat.easy_check$TXHISTORYSTATUS[i] == "Post & Post-adjuvant" & sample.info.dat.easy_check$variant.info[i] == "has variant") {
    sample.info.dat.easy_check$x_numeric[i] <- 9
  }
}

sample.info.dat.easy_check$col <- rep(NA, dim(sample.info.dat)[1])
sample.info.dat.easy_check$col[sample.info.dat.easy_check$variant.info == "freshly frozen"] <- "blue"
sample.info.dat.easy_check$col[sample.info.dat.easy_check$variant.info == "no variant"] <- "red"
sample.info.dat.easy_check$col[sample.info.dat.easy_check$variant.info == "has variant"] <- "black"

sample.info.dat.easy_check$point <- rep(NA, dim(sample.info.dat)[1])
sample.info.dat.easy_check$point[sample.info.dat.easy_check$variant.info == "freshly frozen"] <- 16
sample.info.dat.easy_check$point[sample.info.dat.easy_check$variant.info == "no variant"] <- 17
sample.info.dat.easy_check$point[sample.info.dat.easy_check$variant.info == "has variant"] <- 18

x_num.groups <- sample.info.dat.easy_check %>% group_by(x_numeric) %>% summarise(ct=n())
sample.info.dat.easy_check$x_numbers <- rep(NA, dim(sample.info.dat.easy_check)[1])
for (k in 1:dim(x_num.groups)[1]) {
  cur_num <- as.numeric(x_num.groups[k, 1]); cur_num.ct <- as.numeric(x_num.groups[k, 2])
  if (cur_num.ct == 1) {sample.info.dat.easy_check$x_numbers[sample.info.dat.easy_check$x_numeric == cur_num] <- sample.info.dat.easy_check$x_numeric[sample.info.dat.easy_check$x_numeric == cur_num]}
  else {
    sample.info.dat.easy_check$x_numbers[sample.info.dat.easy_check$x_numeric == cur_num] <- cur_num + seq(-0.25,0.25, by = 0.5/cur_num.ct)[-1]
  }
}

#####################################################################################################################################
#--------------------------- 1. 540 nm ---------------------------
concen.540 <- read.table(file = filepath.concentration.540nm, header = TRUE, colClasses = c("numeric", "numeric"))
colnames(concen.540)[1] <- "BARCODE"
#concen.540 <-  concen.540 %>% drop_na()
dat.540.easy_check <- sample.info.dat.easy_check %>% left_join(concen.540, by="BARCODE")
dat.540.easy_check$remove <- rep(FALSE, dim(dat.540.easy_check)[1])
for (i in 1:dim(dat.540.easy_check)[1]) {
  if (is.na(dat.540.easy_check$conc[i])) {
    dat.540.easy_check$remove[i] = TRUE
  }
}
dat.540.easy_check <- dat.540.easy_check[dat.540.easy_check$remove == FALSE, ]


setwd("/user/xinyugao/ELISA_project/")
pdf("4th_ELISA_450nm-540nm.pdf")
plot(dat.540.easy_check$x_numbers, dat.540.easy_check$conc, xlim = c(0,9.5),ylim = c(-20, max(dat.540.easy_check$conc)+10), 
     ylab="Concentration", xlab="TXHISTORYSTATUS",xaxt = "n", cex=1.2,
     main = "ELISA_450nm-540nm",
     col=dat.540.easy_check$col, pch=dat.540.easy_check$point)
# add text
text(1, -10, "freshly frozen",cex = 0.6, col = "blue"); text(2, -10, "no variant",cex = 0.6, col="red"); text(3, -10, "has variant",cex = 0.6, col="black")
text(5, -10, "no variant",cex = 0.6, col="red");
text(8, -10, "no variant",cex = 0.6, col="red"); #text(9, -10, "has variant",cex = 0.6, col="black")
text(2, -21, "No prior treatment", cex=0.8); text(5,-21, "Post but no-adjuvant", cex=0.8); text(8,-21, "Post & Post-adjuvant", cex=0.8)
# add lines
abline(h=25, col="steelblue", lty=2)
abline(h=6, col="tomato3",  lty=2)

#segments(x0 = 1, y0 = 0, x1 = 5, y1 = 5.15, col = "darkgreen")
patients.variant.info <- dat.540.easy_check %>% group_by(SubjectID) %>% summarise(ct=n())
patients.variant.info <- patients.variant.info[patients.variant.info$ct > 1, ]

filepath.subject.cols <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/plot_results_data/"
load(file=paste0(filepath.subject.cols ,"1st_subject_colors.RData"))
joined.col.info <- patients.variant.info %>% left_join(out.subject.cols, by = "SubjectID")

mycols <- joined.col.info$col
for (j in 1:dim(patients.variant.info)[1]) {
  cur_ID <- as.character(patients.variant.info[j, 1])
  cur.patient.dat <- dat.540.easy_check[dat.540.easy_check$SubjectID == cur_ID, c("x_numbers", "conc")]
  if (dim(cur.patient.dat)[1] < 3) {segments(x0 = as.numeric(cur.patient.dat[1,1]), y0 = as.numeric(cur.patient.dat[1,2]),
                                             x1 = as.numeric(cur.patient.dat[2,1]), y1 = as.numeric(cur.patient.dat[2,2]), col = mycols[j])}
  else {
    for (t in 1:(dim(cur.patient.dat)[1] -1)) {
      segments(x0 = as.numeric(cur.patient.dat[t,1]), y0 = as.numeric(cur.patient.dat[t,2]),
               x1 = as.numeric(cur.patient.dat[t+1,1]), y1 = as.numeric(cur.patient.dat[t+1,2]), col = mycols[j])
    }
  }
}

dev.off()


#--------------------------- 2. 540 nm 30K ---------------------------
concen.540.30K <- read.table(file = filepath.concentration.540nm.30K, header = TRUE, colClasses = c("numeric", "numeric"))
colnames(concen.540.30K)[1] <- "BARCODE"
#concen.540 <-  concen.540 %>% drop_na()
dat.540.30K.easy_check <- sample.info.dat.easy_check %>% left_join(concen.540.30K, by="BARCODE")
dat.540.30K.easy_check$remove <- rep(FALSE, dim(dat.540.30K.easy_check)[1])
for (i in 1:dim(dat.540.30K.easy_check)[1]) {
  if (is.na(dat.540.30K.easy_check$conc[i])) {
    dat.540.30K.easy_check$remove[i] = TRUE
  }
}
dat.540.30K.easy_check <- dat.540.30K.easy_check[dat.540.30K.easy_check$remove == FALSE, ]

setwd("/user/xinyugao/ELISA_project/")
pdf("4th_ELISA_450nm-540nm_30K.pdf")
plot(dat.540.30K.easy_check$x_numbers, dat.540.30K.easy_check$conc, xlim = c(0,9.5),ylim = c(-20, max(dat.540.30K.easy_check$conc)+10), 
     ylab="Concentration", xlab="TXHISTORYSTATUS",xaxt = "n", cex=1.2,
     main = "ELISA_450nm-540nm (30K)",
     col=dat.540.30K.easy_check$col, pch=dat.540.30K.easy_check$point)
# add text
text(1, -10, "freshly frozen",cex = 0.6, col = "blue"); text(2, -10, "no variant",cex = 0.6, col="red"); text(3, -10, "has variant",cex = 0.6, col="black")
text(5, -10, "no variant",cex = 0.6, col="red");
text(8, -10, "no variant",cex = 0.6, col="red"); #text(9, -10, "has variant",cex = 0.6, col="black")
text(2, -21, "No prior treatment", cex=0.8); text(5,-21, "Post but no-adjuvant", cex=0.8); text(8,-21, "Post & Post-adjuvant", cex=0.8)
# add lines
abline(h=25, col="steelblue", lty=2)
abline(h=6, col="tomato3",  lty=2)

#segments(x0 = 1, y0 = 0, x1 = 5, y1 = 5.15, col = "darkgreen")
patients.variant.info <- dat.540.30K.easy_check %>% group_by(SubjectID) %>% summarise(ct=n())
patients.variant.info <- patients.variant.info[patients.variant.info$ct > 1, ]
joined.col.info <- patients.variant.info %>% left_join(out.subject.cols, by = "SubjectID")
mycols <- joined.col.info$col
for (j in 1:dim(patients.variant.info)[1]) {
  cur_ID <- as.character(patients.variant.info[j, 1])
  cur.patient.dat <- dat.540.30K.easy_check[dat.540.30K.easy_check$SubjectID == cur_ID, c("x_numbers", "conc")]
  if (dim(cur.patient.dat)[1] < 3) {segments(x0 = as.numeric(cur.patient.dat[1,1]), y0 = as.numeric(cur.patient.dat[1,2]),
                                             x1 = as.numeric(cur.patient.dat[2,1]), y1 = as.numeric(cur.patient.dat[2,2]), col = mycols[j])}
  else {
    for (t in 1:(dim(cur.patient.dat)[1] -1)) {
      segments(x0 = as.numeric(cur.patient.dat[t,1]), y0 = as.numeric(cur.patient.dat[t,2]),
               x1 = as.numeric(cur.patient.dat[t+1,1]), y1 = as.numeric(cur.patient.dat[t+1,2]), col = mycols[j])
    }
  }
}

dev.off()


#--------------------------- 3. 570 nm ---------------------------
concen.570 <- read.table(file = filepath.concentration.570nm, header = TRUE, colClasses = c("numeric", "numeric"))
colnames(concen.570)[1] <- "BARCODE"
#concen.570 <-  concen.570 %>% drop_na()
dat.570.easy_check <- sample.info.dat.easy_check %>% left_join(concen.570, by="BARCODE")
dat.570.easy_check$remove <- rep(FALSE, dim(dat.570.easy_check)[1])
for (i in 1:dim(dat.570.easy_check)[1]) {
  if (is.na(dat.570.easy_check$conc[i])) {
    dat.570.easy_check$remove[i] = TRUE
  }
}
dat.570.easy_check <- dat.570.easy_check[dat.570.easy_check$remove == FALSE, ]


setwd("/user/xinyugao/ELISA_project/")
pdf("4th_ELISA_450nm-570nm.pdf")
plot(dat.570.easy_check$x_numbers, dat.570.easy_check$conc, xlim = c(0,9.5),ylim = c(-20, max(dat.570.easy_check$conc)+10), 
     ylab="Concentration", xlab="TXHISTORYSTATUS",xaxt = "n", cex=1.2,
     main = "ELISA_450nm-570nm",
     col=dat.570.easy_check$col, pch=dat.570.easy_check$point)
# add text
text(1, -10, "freshly frozen",cex = 0.6, col = "blue"); text(2, -10, "no variant",cex = 0.6, col="red"); text(3, -10, "has variant",cex = 0.6, col="black")
text(5, -10, "no variant",cex = 0.6, col="red");
text(8, -10, "no variant",cex = 0.6, col="red"); #text(9, -10, "has variant",cex = 0.6, col="black")
text(2, -21, "No prior treatment", cex=0.8); text(5,-21, "Post but no-adjuvant", cex=0.8); text(8,-21, "Post & Post-adjuvant", cex=0.8)
# add lines
abline(h=25, col="steelblue", lty=2)
abline(h=6, col="tomato3",  lty=2)

#segments(x0 = 1, y0 = 0, x1 = 5, y1 = 5.15, col = "darkgreen")
patients.variant.info <- dat.570.easy_check %>% group_by(SubjectID) %>% summarise(ct=n())
patients.variant.info <- patients.variant.info[patients.variant.info$ct > 1, ]
joined.col.info <- patients.variant.info %>% left_join(out.subject.cols, by = "SubjectID")
mycols <- joined.col.info$col
for (j in 1:dim(patients.variant.info)[1]) {
  cur_ID <- as.character(patients.variant.info[j, 1])
  cur.patient.dat <- dat.570.easy_check[dat.570.easy_check$SubjectID == cur_ID, c("x_numbers", "conc")]
  if (dim(cur.patient.dat)[1] < 3) {segments(x0 = as.numeric(cur.patient.dat[1,1]), y0 = as.numeric(cur.patient.dat[1,2]),
                                             x1 = as.numeric(cur.patient.dat[2,1]), y1 = as.numeric(cur.patient.dat[2,2]), col = mycols[j])}
  else {
    for (t in 1:(dim(cur.patient.dat)[1] -1)) {
      segments(x0 = as.numeric(cur.patient.dat[t,1]), y0 = as.numeric(cur.patient.dat[t,2]),
               x1 = as.numeric(cur.patient.dat[t+1,1]), y1 = as.numeric(cur.patient.dat[t+1,2]), col = mycols[j])
    }
  }
}

dev.off()


#--------------------------- 4. 570 nm 30K ---------------------------
concen.570.30K <- read.table(file = filepath.concentration.570nm.30K, header = TRUE, colClasses = c("numeric", "numeric"))
colnames(concen.570.30K)[1] <- "BARCODE"
#concen.570 <-  concen.570 %>% drop_na()
dat.570.30K.easy_check <- sample.info.dat.easy_check %>% left_join(concen.570.30K, by="BARCODE")
dat.570.30K.easy_check$remove <- rep(FALSE, dim(dat.570.30K.easy_check)[1])
for (i in 1:dim(dat.570.30K.easy_check)[1]) {
  if (is.na(dat.570.30K.easy_check$conc[i])) {
    dat.570.30K.easy_check$remove[i] = TRUE
  }
}
dat.570.30K.easy_check <- dat.570.30K.easy_check[dat.570.30K.easy_check$remove == FALSE, ]

setwd("/user/xinyugao/ELISA_project/")
pdf("4th_ELISA_450nm-570nm_30K.pdf")
plot(dat.570.30K.easy_check$x_numbers, dat.570.30K.easy_check$conc, xlim = c(0,9.5),ylim = c(-20, max(dat.570.30K.easy_check$conc)+10), 
     ylab="Concentration", xlab="TXHISTORYSTATUS",xaxt = "n", cex=1.2,
     main = "ELISA_450nm-570nm (30K)",
     col=dat.570.30K.easy_check$col, pch=dat.570.30K.easy_check$point)
# add text
text(1, -10, "freshly frozen",cex = 0.6, col = "blue"); text(2, -10, "no variant",cex = 0.6, col="red"); text(3, -10, "has variant",cex = 0.6, col="black")
text(5, -10, "no variant",cex = 0.6, col="red");
text(8, -10, "no variant",cex = 0.6, col="red"); #text(9, -10, "has variant",cex = 0.6, col="black")
text(2, -21, "No prior treatment", cex=0.8); text(5,-21, "Post but no-adjuvant", cex=0.8); text(8,-21, "Post & Post-adjuvant", cex=0.8)
# add lines
abline(h=25, col="steelblue", lty=2)
abline(h=6, col="tomato3",  lty=2)

#segments(x0 = 1, y0 = 0, x1 = 5, y1 = 5.15, col = "darkgreen")
patients.variant.info <- dat.570.30K.easy_check %>% group_by(SubjectID) %>% summarise(ct=n())
patients.variant.info <- patients.variant.info[patients.variant.info$ct > 1, ]
joined.col.info <- patients.variant.info %>% left_join(out.subject.cols, by = "SubjectID")
mycols <- joined.col.info$col
for (j in 1:dim(patients.variant.info)[1]) {
  cur_ID <- as.character(patients.variant.info[j, 1])
  cur.patient.dat <- dat.570.30K.easy_check[dat.570.30K.easy_check$SubjectID == cur_ID, c("x_numbers", "conc")]
  if (dim(cur.patient.dat)[1] < 3) {segments(x0 = as.numeric(cur.patient.dat[1,1]), y0 = as.numeric(cur.patient.dat[1,2]),
                                             x1 = as.numeric(cur.patient.dat[2,1]), y1 = as.numeric(cur.patient.dat[2,2]), col = mycols[j])}
  else {
    for (t in 1:(dim(cur.patient.dat)[1] -1)) {
      segments(x0 = as.numeric(cur.patient.dat[t,1]), y0 = as.numeric(cur.patient.dat[t,2]),
               x1 = as.numeric(cur.patient.dat[t+1,1]), y1 = as.numeric(cur.patient.dat[t+1,2]), col = mycols[j])
    }
  }
}

dev.off()

#-----------------------------------------------------------------------------------------------
# store data
filepath.out <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/plot_results_data/"
out <- list()
out$dat.540 <- dat.540.easy_check
out$dat.540.30K <- dat.540.30K.easy_check
out$dat.570 <- dat.570.easy_check
out$dat.570.30K <- dat.570.30K.easy_check
save(out, file=paste0(filepath.out,"4th_res.RData"))
#-----------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------
# please also create scatter plot between 4thELISA_450nm-570nm and 4thELISA_450nm-570nm.30K, as well as between 4thELISA_450nm-540nm and 4thELISA_450nm-540nm.30K. 
# Please add the regression line and send me the regression results



# 540 and 540 30K
res.4th.540 <- dat.540.30K.easy_check %>% inner_join(dat.540.easy_check, by = c("SubjectID", "x_numeric"))
con.y <- res.4th.540$conc.x;  con.x <- res.4th.540$conc.y

cor(con.x, con.y, method = c("pearson"))
cor(con.x, con.y, method = c("spearman"))

summary(lm(con.y ~ con.x))
setwd("/user/xinyugao/ELISA_project/")
pdf("compare_concentration_4th_540nm_540nm30K.pdf")
plot(con.x, con.y, pch=19, col="black",
     xlab="Concentration (450-540nm)", ylab="Concentration (450-540nm.30K)",
     main = "scatterplot of 4th concentrations (450-540nm and 450-540nm.30K)")
abline(lm(con.y ~ con.x), col = "red")
dev.off()


# 570 and 570 30K
res.4th.570 <- dat.570.30K.easy_check %>% inner_join(dat.570.easy_check, by = c("SubjectID", "x_numeric"))
con.y <- res.4th.570$conc.x;  con.x <- res.4th.570$conc.y
cor(con.x, con.y, method = c("pearson"))
cor(con.x, con.y, method = c("spearman"))


summary(lm(con.y ~ con.x))
setwd("/user/xinyugao/ELISA_project/")
pdf("compare_concentration_4th_570nm_570nm30K.pdf")
plot(con.x, con.y, pch=19, col="black",
     xlab="Concentration (450-570nm)", ylab="Concentration (450-570nm.30K)",
     main = "scatterplot of 4th concentrations (450-570nm and 450-570nm.30K)")
abline(lm(con.y ~ con.x), col = "red")
dev.off()



#-----------------------------------------------------------------------------------------------
# Please also create the scatter plot between the current concentration values and the values for the plasma samples in the 3rd run,  
# add the regression line and send me the regression results
filepath <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/plot_results_data/"
load(file=paste0(filepath, "3rd_res.RData"))
dat.540.plasma.3rd <- out$plasma_540
dat.570.plasma.3rd <- out$plasma_570

res.3rd_4th.540 <- dat.540.plasma.3rd %>% inner_join(dat.540.easy_check, by = c("SubjectID", "x_numeric"))
res.3rd_4th.540.30K <- dat.540.plasma.3rd %>% inner_join(dat.540.30K.easy_check, by = c("SubjectID", "x_numeric"))
res.3rd_4th.570 <- dat.570.plasma.3rd %>% inner_join(dat.570.easy_check, by = c("SubjectID", "x_numeric"))
res.3rd_4th.570.30K <- dat.570.plasma.3rd %>% inner_join(dat.570.30K.easy_check, by = c("SubjectID", "x_numeric"))

# 540
con.y <- res.3rd_4th.540$conc.y;  con.x <- res.3rd_4th.540$conc.x
cor(con.x, con.y, method = c("pearson"))
cor(con.x, con.y, method = c("spearman"))

summary(lm(con.y ~ con.x))
setwd("/user/xinyugao/ELISA_project/")
pdf("compare_concentration_3rd_4th_540nm.pdf")
plot(con.x, con.y, pch=19, col="black",
     xlab="3rd Concentration (450-540nm plasma)", ylab="4th Concentration",
     main = "scatterplot of 3rd (plasma) and 4th concentrations (450-540nm)")
abline(lm(con.y ~ con.x), col = "red")
dev.off()

# 540.30K
con.y <- res.3rd_4th.540.30K$conc.y;  con.x <- res.3rd_4th.540.30K$conc.x
cor(con.x, con.y, method = c("pearson"))
cor(con.x, con.y, method = c("spearman"))

summary(lm(con.y ~ con.x))
setwd("/user/xinyugao/ELISA_project/")
pdf("compare_concentration_3rd_4th_540nm.30K.pdf")
plot(con.x, con.y, pch=19, col="black",
     xlab="3rd Concentration (450-540nm plasma)", ylab="4th Concentration (450-540nm.30K)",
     main = "scatterplot of 3rd (plasma) and 4th (450-540nm.30K) concentrations (450-540nm)")
abline(lm(con.y ~ con.x), col = "red")
dev.off()


# 570
con.y <- res.3rd_4th.570$conc.y;  con.x <- res.3rd_4th.570$conc.x
cor(con.x, con.y, method = c("pearson"))
cor(con.x, con.y, method = c("spearman"))

summary(lm(con.y ~ con.x))
setwd("/user/xinyugao/ELISA_project/")
pdf("compare_concentration_3rd_4th_570nm.pdf")
plot(con.x, con.y, pch=19, col="black",
     xlab="3rd Concentration (450-570nm plasma)", ylab="4th Concentration",
     main = "scatterplot of 3rd (plasma) and 4th concentrations (450-570nm)")
abline(lm(con.y ~ con.x), col = "red")
dev.off()

# 570.30K
con.y <- res.3rd_4th.570.30K$conc.y;  con.x <- res.3rd_4th.570.30K$conc.x
cor(con.x, con.y, method = c("pearson"))
cor(con.x, con.y, method = c("spearman"))

summary(lm(con.y ~ con.x))
setwd("/user/xinyugao/ELISA_project/")
pdf("compare_concentration_3rd_4th_570nm.30K.pdf")
plot(con.x, con.y, pch=19, col="black",
     xlab="3rd Concentration (450-570nm plasma)", ylab="4th Concentration (450-570nm.30K)",
     main = "scatterplot of 3rd (plasma) and 4th (450-570nm.30K) concentrations (450-570nm)")
abline(lm(con.y ~ con.x), col = "red")
dev.off()


#-----------------------------------------------------------------------------------------------
# Please also create the scatter plot between the current concentration values and the values in the 1st run for the serum samples,
# add the regression line and send me the regression results

filepath <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/plot_results_data/"
load(file=paste0(filepath, "1st_res.RData"))
dat.1st <- dat.570
dat.1st.less <- dat.1st[, c("BARCODE", "conc.x", "conc.y", "SubjectID" , "x_numeric")]
dat.1st.540 <- dat.1st.less[, c("BARCODE", "conc.x", "SubjectID" , "x_numeric")]
dat.1st.570 <- dat.1st.less[, c("BARCODE", "conc.y", "SubjectID" , "x_numeric")]

res.1st_4th.540 <- dat.1st.540 %>% inner_join(dat.540.easy_check, by = c("SubjectID", "x_numeric"))
res.1st_4th.540.30K <- dat.1st.540 %>% inner_join(dat.540.30K.easy_check, by = c("SubjectID", "x_numeric"))
res.1st_4th.570 <- dat.1st.570 %>% inner_join(dat.570.easy_check, by = c("SubjectID", "x_numeric"))
res.1st_4th.570.30K <- dat.1st.570 %>% inner_join(dat.570.30K.easy_check, by = c("SubjectID", "x_numeric"))


# 540
con.y <- res.1st_4th.540$conc;  con.x <- res.1st_4th.540$conc.x
cor(con.x, con.y, method = c("pearson"))
cor(con.x, con.y, method = c("spearman"))

summary(lm(con.y ~ con.x))
setwd("/user/xinyugao/ELISA_project/")
pdf("compare_concentration_1st_4th_540nm.pdf")
plot(con.x, con.y, pch=19, col="black",
     xlab="1st Concentration", ylab="4th Concentration",
     main = "scatterplot of 1st and 4th concentrations (450-540nm)")
abline(lm(con.y ~ con.x), col = "red")
dev.off()


# 540.30K
con.y <- res.1st_4th.540.30K$conc;  con.x <- res.1st_4th.540.30K$conc.x
cor(con.x, con.y, method = c("pearson"))
cor(con.x, con.y, method = c("spearman"))

summary(lm(con.y ~ con.x))
setwd("/user/xinyugao/ELISA_project/")
pdf("compare_concentration_1st_4th_540nm.30K.pdf")
plot(con.x, con.y, pch=19, col="black",
     xlab="1st Concentration", ylab="4th Concentration",
     main = "scatterplot of 1st and 4th concentrations (450-540nm.30K)")
abline(lm(con.y ~ con.x), col = "red")
dev.off()


# 570
con.y <- res.1st_4th.570$conc;  con.x <- res.1st_4th.570$conc.y
cor(con.x, con.y, method = c("pearson"))
cor(con.x, con.y, method = c("spearman"))

summary(lm(con.y ~ con.x))
setwd("/user/xinyugao/ELISA_project/")
pdf("compare_concentration_1st_4th_570nm.pdf")
plot(con.x, con.y, pch=19, col="black",
     xlab="1st Concentration", ylab="4th Concentration",
     main = "scatterplot of 1st and 4th concentrations (450-570nm)")
abline(lm(con.y ~ con.x), col = "red")
dev.off()


# 570.30K
con.y <- res.1st_4th.570.30K$conc;  con.x <- res.1st_4th.570.30K$conc.y
cor(con.x, con.y, method = c("pearson"))
cor(con.x, con.y, method = c("spearman"))

summary(lm(con.y ~ con.x))
setwd("/user/xinyugao/ELISA_project/")
pdf("compare_concentration_1st_4th_570nm.30K.pdf")
plot(con.x, con.y, pch=19, col="black",
     xlab="1st Concentration", ylab="4th Concentration",
     main = "scatterplot of 1st and 4th concentrations (450-570nm.30K)")
abline(lm(con.y ~ con.x), col = "red")
dev.off()