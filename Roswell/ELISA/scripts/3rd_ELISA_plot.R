rm(list=ls())
library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr)

# please draw plots for serum and plasma samples respectively using 3rdELISA_450nm-570nm, 3rdELISA_450nm-540nm. 
# Sample information is in "Sheet1" sheet of 3rd_ELISA_samples_5-5-22.xlsx


filepath.concentration.540nm <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/inputs/3rd/3rdELISA_450nm-540nm"
filepath.concentration.570nm <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/inputs/3rd/3rdELISA_450nm-570nm"
filepath.sample.info <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/inputs/3rd/3rd_ELISA_samples_5-5-22.xlsx"

# sample info
sample.info <- read_excel(path = filepath.sample.info, sheet = "Sheet1")
selected.columns <- c("SubjectID", "BARCODE", "TXHISTORYSTATUS" , "group", "rs6494889", "Sample type")
sample.info.dat <- sample.info[, selected.columns]
dim(sample.info.dat)
# 29 * 5
apply(sample.info.dat, 2, function(x) length(unique(x)))
# SubjectID         BARCODE TXHISTORYSTATUS           group       rs6494889 
#     27              27               5               3               4 
table(sample.info.dat$BARCODE)[table(sample.info.dat$BARCODE) > 1]
# 093381175 163371091 
# 2         2

sample.info.dat[sample.info.dat$BARCODE == "093381175", ]
# PT-00216317 093381175 NO PRIOR TREATMENT (EXCLUDING DIAGNOSTIC BIOPSY) Taqman singletimepoint   2 serum        
# PT-00216317 093381175 NO PRIOR TREATMENT (EXCLUDING DIAGNOSTIC BIOPSY) Taqman singletimepoint   2 plasma     
sample.info.dat[sample.info.dat$BARCODE == "163371091", ]
# PT-00342718 163371091 NO PRIOR TREATMENT (EXCLUDING DIAGNOSTIC BIOPSY) Taqman singletimepoint   1 serum        
# PT-00342718 163371091 NO PRIOR TREATMENT (EXCLUDING DIAGNOSTIC BIOPSY) Taqman singletimepoint   1 plasma       

# draw plots for serum and plasma samples respectively
table(sample.info.dat$`Sample type`)
# plasma  serum 
# 10     19 

#---------------------------------------------  pre-works  starts  ----------------------------------------------------------#
table(sample.info.dat$TXHISTORYSTATUS)
# No prior treatment (excluding diagnostic biopsy)    NO PRIOR TREATMENT (EXCLUDING DIAGNOSTIC BIOPSY) POST-SURGICAL AND POST-ADJUVANT OR SYSTEMIC THERAPY 
# 7                                                   8                                                   7 
# Post-surgical/No adjuvant or systemic therapy       POST-SURGICAL/NO ADJUVANT OR SYSTEMIC THERAPY 
# 2                                                   5 

sample.info.dat$TXHISTORYSTATUS[sample.info.dat$TXHISTORYSTATUS == "No prior treatment (excluding diagnostic biopsy)" | 
                                  sample.info.dat$TXHISTORYSTATUS == "NO PRIOR TREATMENT (EXCLUDING DIAGNOSTIC BIOPSY)"] <- "No Prior"
sample.info.dat$TXHISTORYSTATUS[sample.info.dat$TXHISTORYSTATUS == "Post-surgical/No adjuvant or systemic therapy"  |
                                  sample.info.dat$TXHISTORYSTATUS ==  "POST-SURGICAL/NO ADJUVANT OR SYSTEMIC THERAPY"] <- "Post but no-adjuvant"
sample.info.dat$TXHISTORYSTATUS[sample.info.dat$TXHISTORYSTATUS == "POST-SURGICAL AND POST-ADJUVANT OR SYSTEMIC THERAPY"] <- "Post & Post-adjuvant"
table(sample.info.dat$TXHISTORYSTATUS)
# No Prior Post but no-adjuvant Post & Post-adjuvant 
# 15                    7                    7 

sample.info.dat$variant.info <- rep(NA, dim(sample.info.dat)[1])
for (i in 1:dim(sample.info.dat)[1]) {
  if (sample.info.dat$group[i] == "freshly frozen") {sample.info.dat$variant.info[i] <- "freshly frozen"}
  if (!is.na(sample.info.dat$rs6494889[i]) & sample.info.dat$rs6494889[i] == 0) {sample.info.dat$variant.info[i] <- "no variant"}
  if (!is.na(sample.info.dat$rs6494889[i]) & (sample.info.dat$rs6494889[i] == 1 | sample.info.dat$rs6494889[i] == 2)) {sample.info.dat$variant.info[i] <- "has variant"}
}
table(sample.info.dat$variant.info)
# freshly frozen    has variant     no variant 
#   4                    12             13 


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


#---------------------------------------------  pre-works  ends  ----------------------------------------------------------#

sample.info.dat.plasma <- sample.info.dat.easy_check %>% filter(`Sample type` == 'plasma')
dim(sample.info.dat.plasma)
# 10  6
apply(sample.info.dat.plasma, 2, function(x) length(unique(x)))
# SubjectID     BARCODE TXHISTORYSTATUS           group       rs6494889     Sample type 
# 10              10               3               3               4               1
# unique BARCODE

sample.info.dat.serum <- sample.info.dat.easy_check %>% filter(`Sample type` == 'serum')
dim(sample.info.dat.serum)
# 19 6
apply(sample.info.dat.serum, 2, function(x) length(unique(x)))
# SubjectID     BARCODE TXHISTORYSTATUS           group       rs6494889     Sample type 
# 19              19               5               3               4               1
# unique BARCODE
#####################################################################################################################################
#--------------------------- 1. 540 nm ---------------------------
concen.540 <- read.table(file = filepath.concentration.540nm, header = TRUE, colClasses = c("character", "numeric", "character"))
colnames(concen.540)[1] <- "BARCODE"


concen.540.plasma <- concen.540 %>% filter(type == "Plasma")
concen.540.serum <- concen.540 %>% filter(type == "serum")
dat.540.plasma <- sample.info.dat.plasma %>% left_join(concen.540.plasma, by="BARCODE")
dat.540.serum <- sample.info.dat.serum %>% left_join(concen.540.serum, by="BARCODE")



#--------------------------- 1.1 540 nm plasma ---------------------------
x_num.groups <- dat.540.plasma %>% group_by(x_numeric) %>% summarise(ct=n())
dat.540.plasma$x_numbers <- rep(NA, dim(dat.540.plasma)[1])
for (k in 1:dim(x_num.groups)[1]) {
  cur_num <- as.numeric(x_num.groups[k, 1]); cur_num.ct <- as.numeric(x_num.groups[k, 2])
  if (cur_num.ct == 1) {dat.540.plasma$x_numbers[dat.540.plasma$x_numeric == cur_num] <- dat.540.plasma$x_numeric[dat.540.plasma$x_numeric == cur_num]}
  else {
    dat.540.plasma$x_numbers[dat.540.plasma$x_numeric == cur_num] <- cur_num + seq(-0.25,0.25, by = 0.5/cur_num.ct)[-1]
  }
}

setwd("/user/xinyugao/ELISA_project/")
pdf("3rd_ELISA_450nm-540nm_plasma.pdf")
plot(dat.540.plasma$x_numbers, dat.540.plasma$conc, xlim = c(0,9.5),ylim = c(-20, max(dat.540.plasma$conc)+10), 
     ylab="Concentration", xlab="TXHISTORYSTATUS",xaxt = "n", cex=1.2,
     main = "ELISA_450nm-540nm (plasma)",
     col=dat.540.plasma$col, pch=dat.540.plasma$point)
# add text
text(1, -10, "freshly frozen",cex = 0.6, col = "blue"); text(2, -10, "no variant",cex = 0.6, col="red"); text(3, -10, "has variant",cex = 0.6, col="black")
text(5, -10, "no variant",cex = 0.6, col="red"); 
#text(7, -10, "freshly frozen",cex = 0.6, col = "blue"); text(8, -10, "no variant",cex = 0.6, col="red"); text(9, -10, "has variant",cex = 0.6, col="black")
text(2, -21, "No prior treatment", cex=0.8); text(5,-21, "Post but no-adjuvant", cex=0.8); text(8,-21, "Post & Post-adjuvant", cex=0.8)
# add lines
abline(h=25, col="steelblue", lty=2)
abline(h=6, col="tomato3",  lty=2)
#segments(x0 = 1, y0 = 0, x1 = 5, y1 = 5.15, col = "darkgreen")
dev.off()




#--------------------------- 1.2 540 serum ---------------------------
x_num.groups <- dat.540.serum %>% group_by(x_numeric) %>% summarise(ct=n())
dat.540.serum$x_numbers <- rep(NA, dim(dat.540.serum)[1])
for (k in 1:dim(x_num.groups)[1]) {
  cur_num <- as.numeric(x_num.groups[k, 1]); cur_num.ct <- as.numeric(x_num.groups[k, 2])
  if (cur_num.ct == 1) {dat.540.serum$x_numbers[dat.540.serum$x_numeric == cur_num] <- dat.540.serum$x_numeric[dat.540.serum$x_numeric == cur_num]}
  else {
    dat.540.serum$x_numbers[dat.540.serum$x_numeric == cur_num] <- cur_num + seq(-0.25,0.25, by = 0.5/cur_num.ct)[-1]
  }
}


setwd("/user/xinyugao/ELISA_project/")
pdf("3rd_ELISA_450nm-540nm_serum.pdf")
plot(dat.540.serum$x_numbers, dat.540.serum$conc, xlim = c(0,9.5),ylim = c(-20, max(dat.540.serum$conc)+10), 
     ylab="Concentration", xlab="TXHISTORYSTATUS",xaxt = "n", cex=1.2,
     main = "ELISA_450nm-540nm (serum)",
     col=dat.540.serum$col, pch=dat.540.serum$point)
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



#--------------------------- 2. 570 nm ---------------------------
concen.570 <- read.table(file = filepath.concentration.570nm, header = TRUE, colClasses = c("character", "numeric", "character"))
colnames(concen.570)[1] <- "BARCODE"

concen.570.plasma <- concen.570 %>% filter(type == "Plasma")
concen.570.serum <- concen.570 %>% filter(type == "serum")
dat.570.plasma <- sample.info.dat.plasma %>% left_join(concen.570.plasma, by="BARCODE")
dat.570.serum <- sample.info.dat.serum %>% left_join(concen.570.serum, by="BARCODE")

#--------------------------- 2.1 570 nm plasma ---------------------------
x_num.groups <- dat.570.plasma %>% group_by(x_numeric) %>% summarise(ct=n())
dat.570.plasma$x_numbers <- rep(NA, dim(dat.570.plasma)[1])
for (k in 1:dim(x_num.groups)[1]) {
  cur_num <- as.numeric(x_num.groups[k, 1]); cur_num.ct <- as.numeric(x_num.groups[k, 2])
  if (cur_num.ct == 1) {dat.570.plasma$x_numbers[dat.570.plasma$x_numeric == cur_num] <- dat.570.plasma$x_numeric[dat.570.plasma$x_numeric == cur_num]}
  else {
    dat.570.plasma$x_numbers[dat.570.plasma$x_numeric == cur_num] <- cur_num + seq(-0.25,0.25, by = 0.5/cur_num.ct)[-1]
  }
}

setwd("/user/xinyugao/ELISA_project/")
pdf("3rd_ELISA_450nm-570nm_plasma.pdf")
plot(dat.570.plasma$x_numbers, dat.570.plasma$conc, xlim = c(0,9.5),ylim = c(-20, max(dat.570.plasma$conc)+10), 
     ylab="Concentration", xlab="TXHISTORYSTATUS",xaxt = "n", cex=1.2,
     main = "ELISA_450nm-570nm (plasma)",
     col=dat.570.plasma$col, pch=dat.570.plasma$point)
# add text
text(1, -10, "freshly frozen",cex = 0.6, col = "blue"); text(2, -10, "no variant",cex = 0.6, col="red"); text(3, -10, "has variant",cex = 0.6, col="black")
text(5, -10, "no variant",cex = 0.6, col="red"); 
#text(7, -10, "freshly frozen",cex = 0.6, col = "blue"); text(8, -10, "no variant",cex = 0.6, col="red"); text(9, -10, "has variant",cex = 0.6, col="black")
text(2, -21, "No prior treatment", cex=0.8); text(5,-21, "Post but no-adjuvant", cex=0.8); text(8,-21, "Post & Post-adjuvant", cex=0.8)
# add lines
abline(h=25, col="steelblue", lty=2)
abline(h=6, col="tomato3",  lty=2)
#segments(x0 = 1, y0 = 0, x1 = 5, y1 = 5.15, col = "darkgreen")
dev.off()



#--------------------------- 2.2 570 serum ---------------------------

x_num.groups <- dat.570.serum %>% group_by(x_numeric) %>% summarise(ct=n())
dat.570.serum$x_numbers <- rep(NA, dim(dat.570.serum)[1])
for (k in 1:dim(x_num.groups)[1]) {
  cur_num <- as.numeric(x_num.groups[k, 1]); cur_num.ct <- as.numeric(x_num.groups[k, 2])
  if (cur_num.ct == 1) {dat.570.serum$x_numbers[dat.570.serum$x_numeric == cur_num] <- dat.570.serum$x_numeric[dat.570.serum$x_numeric == cur_num]}
  else {
    dat.570.serum$x_numbers[dat.570.serum$x_numeric == cur_num] <- cur_num + seq(-0.25,0.25, by = 0.5/cur_num.ct)[-1]
  }
}


setwd("/user/xinyugao/ELISA_project/")
pdf("3rd_ELISA_450nm-570nm_serum.pdf")
plot(dat.570.serum$x_numbers, dat.570.serum$conc, xlim = c(0,9.5),ylim = c(-20, max(dat.570.serum$conc)+10), 
     ylab="Concentration", xlab="TXHISTORYSTATUS",xaxt = "n", cex=1.2,
     main = "ELISA_450nm-570nm (serum)",
     col=dat.570.serum$col, pch=dat.570.serum$point)
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

#--------------------------------------------------------------------------------------------------------------------
# Store results 
filepath.out <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/plot_results_data/"
out <- list()
out$plasma_540 <- dat.540.plasma
out$serum_540 <- dat.540.serum
out$plasma_570 <- dat.570.plasma
out$serum_570 <- dat.570.serum
save(out, file=paste0(filepath.out,"3rd_res.RData"))

############################################################################################################################################
# Task 2: Please also create the scatter plot between the current concentration values and the values in the 1st run for the serum samples, 
#         add the regression line and send me the regression results
filepath <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/plot_results_data/"
load(file=paste0(filepath, "3rd_res.RData"))
dat.540.serum <- out$serum_540
dat.570.serum <- out$serum_570

filepath <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/plot_results_data/"
load(file=paste0(filepath, "1st_res.RData"))
dat.1st <- dat.570
dat.1st.less <- dat.1st[, c("BARCODE", "conc.x", "conc.y", "SubjectID" , "x_numeric")]

# 540
res.1st_3rd.540 <- dat.540.serum %>% inner_join(dat.1st.less, by = c("SubjectID", "x_numeric"))
summary(lm(res.1st_3rd.540$conc ~ res.1st_3rd.540$conc.x))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
#   (Intercept)            111.5163    19.9931   5.578   0.0307 *
#   res.1st_3rd.540$conc.x   0.8542     0.3866   2.209   0.1578  
# Residual standard error: 29.96 on 2 degrees of freedom
# Multiple R-squared:  0.7094,	Adjusted R-squared:  0.564 
# F-statistic: 4.881 on 1 and 2 DF,  p-value: 0.1578

setwd("/user/xinyugao/ELISA_project/")
pdf("compare_concentration_1st-3rd_450-540_serum.pdf")
plot(res.1st_3rd.540$conc.x, res.1st_3rd.540$conc, pch=19, col="black",
     xlab="1st Concentration", ylab="3rd Concentration",
     main = "scatterplot of 3rd and 1st concentrations (450-540nm serum)")
abline(lm(res.1st_3rd.540$conc ~ res.1st_3rd.540$conc.x), col = "red")
dev.off()

# 570
res.1st_3rd.570 <- dat.570.serum %>% inner_join(dat.1st.less, by = c("SubjectID", "x_numeric"))
summary(lm(res.1st_3rd.570$conc ~ res.1st_3rd.570$conc.y))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)            111.9095    20.0193   5.590   0.0305 *
#   res.1st_3rd.570$conc.y   0.8541     0.3876   2.204   0.1584  
# Residual standard error: 30.04 on 2 degrees of freedom
# Multiple R-squared:  0.7083,	Adjusted R-squared:  0.5624 
# F-statistic: 4.856 on 1 and 2 DF,  p-value: 0.1584

setwd("/user/xinyugao/ELISA_project/")
pdf("compare_concentration_1st-3rd_450-570_serum.pdf")
plot(res.1st_3rd.570$conc.y, res.1st_3rd.570$conc, pch=19, col="black",
     xlab="1st Concentration", ylab="3rd Concentration",
     main = "scatterplot of 3rd and 1st concentrations (450-570nm serum)")
abline(lm(res.1st_3rd.570$conc ~ res.1st_3rd.570$conc.y), col = "red")
dev.off()



############################################################################################################################################
# Task 3: Please also create the scatter plot between the current concentration values of the plasma samples 
#         and the concentration values of their corresponding serum samples in the 1st run, 
#         add the regression line and send me the regression results


dat.540.plasma <- out$plasma_540
dat.570.plasma <- out$plasma_570

# 540
res.1st_3rd.540 <- dat.540.plasma %>% inner_join(dat.1st.less, by = c("SubjectID", "x_numeric"))
summary(lm(res.1st_3rd.540$conc ~ res.1st_3rd.540$conc.x))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
#   (Intercept)             49.2161    10.2714   4.792  0.00492 **
#   res.1st_3rd.540$conc.x  -0.6912     0.8057  -0.858  0.43017   
# Residual standard error: 22.79 on 5 degrees of freedom
# Multiple R-squared:  0.1283,	Adjusted R-squared:  -0.04604 
# F-statistic: 0.7359 on 1 and 5 DF,  p-value: 0.4302
setwd("/user/xinyugao/ELISA_project/")
pdf("compare_concentration_1st(serum)-3rd(plasma)_450-540.pdf")
plot(res.1st_3rd.540$conc.x, res.1st_3rd.540$conc, pch=19, col="black",
     xlab="1st Concentration", ylab="3rd Concentration",
     main = "scatterplot of 3rd (plasma) and 1st (serum) concentrations (450-540nm)")
abline(lm(res.1st_3rd.540$conc ~ res.1st_3rd.540$conc.x), col = "red")
dev.off()


# 570
res.1st_3rd.570 <- dat.570.plasma %>% inner_join(dat.1st.less, by = c("SubjectID", "x_numeric"))
summary(lm(res.1st_3rd.570$conc ~ res.1st_3rd.570$conc.y))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)             49.2703    10.3310   4.769  0.00502 **
#   res.1st_3rd.570$conc.y  -0.7003     0.8171  -0.857  0.43057   
# Residual standard error: 22.92 on 5 degrees of freedom
# Multiple R-squared:  0.1281,	Adjusted R-squared:  -0.04629 
# F-statistic: 0.7345 on 1 and 5 DF,  p-value: 0.4306
setwd("/user/xinyugao/ELISA_project/")
pdf("compare_concentration_1st(serum)-3rd(plasma)_450-570.pdf")
plot(res.1st_3rd.540$conc.y, res.1st_3rd.540$conc, pch=19, col="black",
     xlab="1st Concentration", ylab="3rd Concentration",
     main = "scatterplot of 3rd (plasma) and 1st (serum) concentrations (450-570nm)")
abline(lm(res.1st_3rd.570$conc ~ res.1st_3rd.570$conc.y), col = "red")
dev.off()




############################################################################################################################################
# Task 4: Please also create the scatter plot between the current concentration values (plasma) and the values in the 2nd run for the serum samples,
#         add the regression line and send me the regression results

filepath <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/plot_results_data/"
load(file=paste0(filepath, "2nd_res.RData"))
dat.2nd <- dat.570
dat.2nd.less <- dat.2nd[, c("BARCODE", "conc.x", "conc.y", "SubjectID" , "x_numeric")]

# 540
res.2nd_3rd.540 <- dat.540.plasma %>% inner_join(dat.2nd.less, by = c("SubjectID", "x_numeric"))
summary(lm(res.2nd_3rd.540$conc ~ res.2nd_3rd.540$conc.x))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)            52.63571   38.38263   1.371    0.401
# res.2nd_3rd.540$conc.x  0.08274    0.59631   0.139    0.912
# Residual standard error: 40.7 on 1 degrees of freedom
# Multiple R-squared:  0.01889,	Adjusted R-squared:  -0.9622 
# F-statistic: 0.01925 on 1 and 1 DF,  p-value: 0.9122
setwd("/user/xinyugao/ELISA_project/")
pdf("compare_concentration_2nd(serum)-3rd(plasma)_450-540.pdf")
plot(res.2nd_3rd.540$conc.x, res.2nd_3rd.540$conc, pch=19, col="black",
     xlab="2nd Concentration", ylab="3rd Concentration",
     main = "scatterplot of 3rd (plasma) and 2nd (serum) concentrations (450-540nm)")
abline(lm(res.2nd_3rd.540$conc ~ res.2nd_3rd.540$conc.x), col = "red")
dev.off()

# 570
res.2nd_3rd.570 <- dat.570.plasma %>% inner_join(dat.2nd.less, by = c("SubjectID", "x_numeric"))
summary(lm(res.2nd_3rd.570$conc ~ res.2nd_3rd.570$conc.y))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)            52.90793   38.41367   1.377    0.400
# res.2nd_3rd.570$conc.y  0.08101    0.59771   0.136    0.914
# Residual standard error: 40.77 on 1 degrees of freedom
# Multiple R-squared:  0.01804,	Adjusted R-squared:  -0.9639 
# F-statistic: 0.01837 on 1 and 1 DF,  p-value: 0.9142
setwd("/user/xinyugao/ELISA_project/")
pdf("compare_concentration_2nd(serum)-3rd(plasma)_450-570.pdf")
plot(res.2nd_3rd.570$conc.y, res.2nd_3rd.570$conc, pch=19, col="black",
     xlab="2nd Concentration", ylab="3rd Concentration",
     main = "scatterplot of 3rd (plasma) and 2nd (serum) concentrations (450-570nm)")
abline(lm(res.2nd_3rd.570$conc ~ res.2nd_3rd.570$conc.y), col = "red")
dev.off()


############################################################################################################################################
# Task 4: Please also create the scatter plot between the current concentration values (serum) and the values in the 2nd run for the serum samples,
#         add the regression line and send me the regression results
# 540
res.2nd_3rd.540 <- dat.540.serum %>% inner_join(dat.2nd.less, by = c("SubjectID", "x_numeric"))
summary(lm(res.2nd_3rd.540$conc ~ res.2nd_3rd.540$conc.x))
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            124.7704    14.2772   8.739 1.07e-07 ***
#   res.2nd_3rd.540$conc.x   0.2267     0.1727   1.313    0.207    
# Residual standard error: 56.81 on 17 degrees of freedom
# Multiple R-squared:  0.09206,	Adjusted R-squared:  0.03865 
# F-statistic: 1.724 on 1 and 17 DF,  p-value: 0.2067

setwd("/user/xinyugao/ELISA_project/")
pdf("compare_concentration_2nd(serum)-3rd(serum)_450-540.pdf")
plot(res.2nd_3rd.540$conc.x, res.2nd_3rd.540$conc, pch=19, col="black",
     xlab="2nd Concentration", ylab="3rd Concentration",
     main = "scatterplot of 3rd (serum) and 2nd (serum) concentrations (450-540nm)")
abline(lm(res.2nd_3rd.540$conc ~ res.2nd_3rd.540$conc.x), col = "red")
dev.off()

# 570
res.2nd_3rd.570 <- dat.570.serum %>% inner_join(dat.2nd.less, by = c("SubjectID", "x_numeric"))
summary(lm(res.2nd_3rd.570$conc ~ res.2nd_3rd.570$conc.y))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
#    (Intercept)            125.0002    14.3015   8.740 1.07e-07 ***
#   res.2nd_3rd.570$conc.y   0.2272     0.1729   1.314    0.206    
# Residual standard error: 56.88 on 17 degrees of freedom
# Multiple R-squared:  0.09224,	Adjusted R-squared:  0.03885 
# F-statistic: 1.728 on 1 and 17 DF,  p-value: 0.2062
setwd("/user/xinyugao/ELISA_project/")
pdf("compare_concentration_2nd(serum)-3rd(serum)_450-570.pdf")
plot(res.2nd_3rd.570$conc.y, res.2nd_3rd.570$conc, pch=19, col="black",
     xlab="2nd Concentration", ylab="3rd Concentration",
     main = "scatterplot of 3rd (serum) and 2nd (serum) concentrations (450-570nm)")
abline(lm(res.2nd_3rd.570$conc ~ res.2nd_3rd.570$conc.y), col = "red")
dev.off()



