rm(list=ls())
library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr)


#------------------------ load 3rd results ------------------------
filepath <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/plot_results_data/"
load(file=paste0(filepath, "3rd_res.RData"))
dat.540.serum <- out$serum_540
dat.570.serum <- out$serum_570
dat.540.plasma <- out$plasma_540
dat.570.plasma <- out$plasma_570


#------------------------ load 1st results ------------------------
filepath <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/plot_results_data/"
load(file=paste0(filepath, "1st_res.RData"))
dat.1st <- dat.570
dat.1st.less <- dat.1st[, c("BARCODE", "conc.x", "conc.y", "SubjectID" , "x_numeric")]

#------------------------ load 2nd results ------------------------
filepath <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/plot_results_data/"
load(file=paste0(filepath, "2nd_res.RData"))
dat.2nd <- dat.570
dat.2nd.less <- dat.2nd[, c("BARCODE", "conc.x", "conc.y", "SubjectID" , "x_numeric")]


#----------------------------------------------------------------------------
# Task 1: calculate the pearson and spearman correlation coefficients between 2nd and 1st run
res.1st.2nd <- dat.1st.less %>% inner_join(dat.2nd.less, by = c("SubjectID", "x_numeric"))
colnames(res.1st.2nd) <- c("Barcode.1st", "conc.1st.570", "conc.1st.540", "SubjectID", "specific_groups", 
                           "Barcode.2nd", "conc.2nd.570", "conc.2nd.540")
# 540
x = res.1st.2nd$conc.1st.540; y = res.1st.2nd$conc.2nd.540
# 7 
cor(x, y, method = c("pearson"))
# 0.976637
cor(x, y, method = c("spearman"))
# 0.8669214

# 570
x = res.1st.2nd$conc.1st.570; y = res.1st.2nd$conc.2nd.570
# 7 
cor(x, y, method = c("pearson"))
# 0.9764312
cor(x, y, method = c("spearman"))
# 0.8669214


#----------------------------------------------------------------------------
# Task 2: for serum samples and plasma samples respectively between 3rd run and 1st run

######  serum
# 540
res.1st_3rd.serum.540 <- dat.540.serum %>% inner_join(dat.1st.less, by = c("SubjectID", "x_numeric"))
x = res.1st_3rd.serum.540$conc; y = res.1st_3rd.serum.540$conc.x
# 4
cor(x, y, method = c("pearson"))
#  0.8422339
cor(x, y, method = c("spearman"))
# 0.8


# 570
res.1st_3rd.serum.570 <- dat.570.serum %>% inner_join(dat.1st.less, by = c("SubjectID", "x_numeric"))
x = res.1st_3rd.serum.570$conc; y = res.1st_3rd.serum.570$conc.y
# 4
cor(x, y, method = c("pearson"))
#  0.841596
cor(x, y, method = c("spearman"))
# 0.8


## plasma
# 540
res.1st_3rd.plasma.540 <- dat.540.plasma %>% inner_join(dat.1st.less, by = c("SubjectID", "x_numeric"))
x = res.1st_3rd.plasma.540$conc; y = res.1st_3rd.plasma.540$conc.x
# 7
cor(x, y, method = c("pearson"))
# -0.3581888
cor(x, y, method = c("spearman"))
# -0.1970276

# 570
res.1st_3rd.plasma.570 <- dat.570.plasma %>% inner_join(dat.1st.less, by = c("SubjectID", "x_numeric"))
x = res.1st_3rd.plasma.570$conc; y = res.1st_3rd.plasma.570$conc.y
# 7
cor(x, y, method = c("pearson"))
#  -0.3578952
cor(x, y, method = c("spearman"))
# -0.1970276



#----------------------------------------------------------------------------
# Task 3: only serum samples between the 3rd run and 2nd run

# 540
res.2nd_3rd.serum.540 <- dat.540.serum %>% inner_join(dat.2nd.less, by = c("SubjectID", "x_numeric"))
# 19
x = res.2nd_3rd.serum.540$conc; y = res.2nd_3rd.serum.540$conc.x
cor(x, y, method = c("pearson"))
#  0.3034168
cor(x, y, method = c("spearman"))
# 0.6440891


# 570
res.2nd_3rd.serum.570 <- dat.570.serum %>% inner_join(dat.2nd.less, by = c("SubjectID", "x_numeric"))
x = res.2nd_3rd.serum.570$conc; y = res.2nd_3rd.serum.570$conc.y
cor(x, y, method = c("pearson"))
#  0.3037184
cor(x, y, method = c("spearman"))
#0.6335866


