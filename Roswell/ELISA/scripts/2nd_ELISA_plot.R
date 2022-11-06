# Could you please draw the plots again using this new data?
# The excel file contains sample information. Please use the "2nd_ELISA_samples" sheet.

rm(list=ls())
library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr)


filepath.concentration.540nm <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/inputs/2ndELISA_450nm-540nm"
filepath.concentration.570nm <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/inputs/2ndELISA_450nm-570nm"
filepath.sample.info <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/inputs/2nd_ELISA_samples_results.xlsx"

# sample info
sample.info <- read_excel(path = filepath.sample.info, sheet = "2nd_ELISA_samples")
selected.columns <- c("SubjectID", "BARCODE", "TXHISTORYSTATUS" , "group", "rs6494889" )
sample.info.dat <- sample.info[, selected.columns]
dim(sample.info.dat)
# 29 * 5
apply(sample.info.dat, 2, function(x) length(unique(x)))
# SubjectID         BARCODE TXHISTORYSTATUS           group       rs6494889 
#     29              29               5               3               4 
# subject ID and barcode are unique columns
# subject ID unqiue -> no need to connect any points

###############################################################################################################################
# 540 nm
concen.540 <- read.table(file = filepath.concentration.540nm, header = TRUE, colClasses = c("character", "numeric"))
colnames(concen.540)[1] <- "BARCODE"
dat.540 <- merge(concen.540, sample.info.dat, by="BARCODE")

table(dat.540$TXHISTORYSTATUS)
# No prior treatment (excluding diagnostic biopsy)    NO PRIOR TREATMENT (EXCLUDING DIAGNOSTIC BIOPSY)    POST-SURGICAL AND POST-ADJUVANT OR SYSTEMIC THERAPY 
#                 4                                                   8                                                   8 
# Post-surgical/No adjuvant or systemic therapy       POST-SURGICAL/NO ADJUVANT OR SYSTEMIC THERAPY 
#                1                                                   8 
dat.540$TXHISTORYSTATUS[dat.540$TXHISTORYSTATUS == "No prior treatment (excluding diagnostic biopsy)" | 
                          dat.540$TXHISTORYSTATUS == "NO PRIOR TREATMENT (EXCLUDING DIAGNOSTIC BIOPSY)"] <- "No Prior"
dat.540$TXHISTORYSTATUS[dat.540$TXHISTORYSTATUS == "Post-surgical/No adjuvant or systemic therapy"  |
                          dat.540$TXHISTORYSTATUS ==  "POST-SURGICAL/NO ADJUVANT OR SYSTEMIC THERAPY"] <- "Post but no-adjuvant"
dat.540$TXHISTORYSTATUS[dat.540$TXHISTORYSTATUS == "POST-SURGICAL AND POST-ADJUVANT OR SYSTEMIC THERAPY"] <- "Post & Post-adjuvant"
table(dat.540$TXHISTORYSTATUS)
# No Prior       Post but no-adjuvant      Post & Post-adjuvant 
#   12                    9                    8 

dat.540$variant.info <- rep(NA, dim(dat.540)[1])
for (i in 1:dim(dat.540)[1]) {
  if (dat.540$group[i] == "freshly frozen") {dat.540$variant.info[i] <- "freshly frozen"}
  if (!is.na(dat.540$rs6494889[i]) & dat.540$rs6494889[i] == 0) {dat.540$variant.info[i] <- "no variant"}
  if (!is.na(dat.540$rs6494889[i]) & (dat.540$rs6494889[i] == 1 | dat.540$rs6494889[i] == 2)) {dat.540$variant.info[i] <- "has variant"}
}

dat.540$TXHISTORYSTATUS <- factor(dat.540$TXHISTORYSTATUS, levels = c("No Prior", "Post but no-adjuvant", 
                                                                      "Post & Post-adjuvant"))
dat.540$variant.info <- factor(dat.540$variant.info, levels = c("freshly frozen", "no variant", "has variant"))
dat.540.easy_check <- dat.540 %>% arrange(TXHISTORYSTATUS, variant.info)



dat.540.easy_check$x_numeric <- rep(NA, dim(dat.540)[1])
for (i in 1:dim(dat.540)[1]) {
  if (dat.540.easy_check$TXHISTORYSTATUS[i] == "No Prior" & dat.540.easy_check$variant.info[i] == "freshly frozen") {
    dat.540.easy_check$x_numeric[i] <- 1
  }
  else if (dat.540.easy_check$TXHISTORYSTATUS[i] == "No Prior" & dat.540.easy_check$variant.info[i] == "no variant") {
    dat.540.easy_check$x_numeric[i] <- 2
  }
  else if (dat.540.easy_check$TXHISTORYSTATUS[i] == "No Prior" & dat.540.easy_check$variant.info[i] == "has variant") {
    dat.540.easy_check$x_numeric[i] <- 3
  }
  else if (dat.540.easy_check$TXHISTORYSTATUS[i] == "Post but no-adjuvant" & dat.540.easy_check$variant.info[i] == "freshly frozen") {
    dat.540.easy_check$x_numeric[i] <- 4
  }
  else if (dat.540.easy_check$TXHISTORYSTATUS[i] == "Post but no-adjuvant" & dat.540.easy_check$variant.info[i] == "no variant") {
    dat.540.easy_check$x_numeric[i] <- 5
  }
  else if (dat.540.easy_check$TXHISTORYSTATUS[i] == "Post but no-adjuvant" & dat.540.easy_check$variant.info[i] == "has variant") {
    dat.540.easy_check$x_numeric[i] <- 6
  }
  else if (dat.540.easy_check$TXHISTORYSTATUS[i] == "Post & Post-adjuvant" & dat.540.easy_check$variant.info[i] == "freshly frozen") {
    dat.540.easy_check$x_numeric[i] <- 7
  }
  else if (dat.540.easy_check$TXHISTORYSTATUS[i] == "Post & Post-adjuvant" & dat.540.easy_check$variant.info[i] == "no variant") {
    dat.540.easy_check$x_numeric[i] <- 8
  }
  else if (dat.540.easy_check$TXHISTORYSTATUS[i] == "Post & Post-adjuvant" & dat.540.easy_check$variant.info[i] == "has variant") {
    dat.540.easy_check$x_numeric[i] <- 9
  }
}

dat.540.easy_check$col <- rep(NA, dim(dat.540)[1])
dat.540.easy_check$col[dat.540.easy_check$variant.info == "freshly frozen"] <- "blue"
dat.540.easy_check$col[dat.540.easy_check$variant.info == "no variant"] <- "red"
dat.540.easy_check$col[dat.540.easy_check$variant.info == "has variant"] <- "black"

dat.540.easy_check$point <- rep(NA, dim(dat.540)[1])
dat.540.easy_check$point[dat.540.easy_check$variant.info == "freshly frozen"] <- 16
dat.540.easy_check$point[dat.540.easy_check$variant.info == "no variant"] <- 17
dat.540.easy_check$point[dat.540.easy_check$variant.info == "has variant"] <- 18

x_num.groups <- dat.540.easy_check %>% group_by(x_numeric) %>% summarise(ct=n())
dat.540.easy_check$x_numbers <- rep(NA, dim(dat.540)[1])
for (k in 1:dim(x_num.groups)[1]) {
  cur_num <- as.numeric(x_num.groups[k, 1]); cur_num.ct <- as.numeric(x_num.groups[k, 2])
  if (cur_num.ct == 1) {dat.540.easy_check$x_numbers[dat.540.easy_check$x_numeric == cur_num] <- dat.540.easy_check$x_numeric[dat.540.easy_check$x_numeric == cur_num]}
  else {
    dat.540.easy_check$x_numbers[dat.540.easy_check$x_numeric == cur_num] <- cur_num + seq(-0.25,0.25, by = 0.5/cur_num.ct)[-1]
  }
}


setwd("/user/xinyugao/ELISA_project/")
pdf("2nd_ELISA_450nm-540nm.pdf")

plot(dat.540.easy_check$x_numbers, dat.540.easy_check$conc, xlim = c(0,9.5),ylim = c(-20, max(dat.540.easy_check$conc)+10), 
     ylab="Concentration", xlab="TXHISTORYSTATUS",xaxt = "n", cex=1.2,
     main = "ELISA_450nm-540nm",
     col=dat.540.easy_check$col, pch=dat.540.easy_check$point)
# add text
text(1, -10, "freshly frozen",cex = 0.6, col = "blue"); text(2, -10, "no variant",cex = 0.6, col="red"); text(3, -10, "has variant",cex = 0.6, col="black")
text(5, -10, "no variant",cex = 0.6, col="red"); text(6, -10, "has variant",cex = 0.6, col="black")
text(8, -10, "no variant",cex = 0.6, col="red"); text(9, -10, "has variant",cex = 0.6, col="black")
text(2, -21, "No prior treatment", cex=0.8); text(5,-21, "Post but no-adjuvant", cex=0.8); text(8,-21, "Post & Post-adjuvant", cex=0.8)
# add lines
abline(h=25, col="steelblue", lty=2)
abline(h=6, col="tomato3",  lty=2)
#segments(x0 = 1, y0 = 0, x1 = 5, y1 = 5.15, col = "darkgreen")

dev.off()






###############################################################################################################################
# 570 nm
concen.570 <- read.table(file = filepath.concentration.570nm, header = TRUE, colClasses = c("character", "numeric"))
colnames(concen.570)[1] <- "BARCODE"
dat.570 <- merge(concen.570, dat.540.easy_check, by="BARCODE")


setwd("/user/xinyugao/ELISA_project/")
pdf("2nd_ELISA_450nm-570nm.pdf")

plot(dat.570$x_numbers, dat.570$conc.x, xlim = c(0,9.5),ylim = c(-20, max(dat.570$conc.x)+10), 
     ylab="Concentration", xlab="TXHISTORYSTATUS",xaxt = "n", cex=1.2,
     main = "ELISA_450nm-570nm",
     col=dat.570$col, pch=dat.570$point)
# add text
text(1, -10, "freshly frozen",cex = 0.6, col = "blue"); text(2, -10, "no variant",cex = 0.6, col="red"); text(3, -10, "has variant",cex = 0.6, col="black")
text(5, -10, "no variant",cex = 0.6, col="red"); text(6, -10, "has variant",cex = 0.6, col="black")
text(8, -10, "no variant",cex = 0.6, col="red"); text(9, -10, "has variant",cex = 0.6, col="black")
text(2, -21, "No prior treatment", cex=0.8); text(5,-21, "Post but no-adjuvant", cex=0.8); text(8,-21, "Post & Post-adjuvant", cex=0.8)
# add lines
abline(h=25, col="steelblue", lty=2)
abline(h=6, col="tomato3",  lty=2)

dev.off()


# store data
filepath.out <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/plot_results_data/"
save(dat.570, file=paste0(filepath.out,"2nd_res.RData"))


