

rm(list=ls())
library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr)





filepath <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/outputs/plot_results_data/"
load(file=paste0(filepath, "1st_res.RData"))
dat.1st <- dat.570

load(file=paste0(filepath, "2nd_res.RData"))
dat.2nd <- dat.570



head(dat.1st)
dat.1st.less <- dat.1st[, c("BARCODE", "conc.x", "conc.y", "SubjectID" , "x_numeric")]
dat.2nd.less <- dat.2nd[, c("BARCODE", "conc.x", "conc.y", "SubjectID" , "x_numeric")]


res.1st.2nd <- dat.1st.less %>% inner_join(dat.2nd.less, by = c("SubjectID", "x_numeric"))
# x_numeric
# no prior treatment:        1                   2               3
#                       freshly frozen     no variant      has variant
# post but no adjuvant:      4                   5               6
#                       freshly frozen     no variant      has variant
# post & post adjuvant:      7                   8               9
#                       freshly frozen     no variant      has variant

colnames(res.1st.2nd) <- c("Barcode.1st", "conc.1st.570", "conc.1st.540", "SubjectID", "specific_groups", 
                           "Barcode.2nd", "conc.2nd.570", "conc.2nd.540")
res.1st.2nd

# scatterplot of 540
# linear regression
summary(lm(res.1st.2nd$conc.2nd.540 ~ res.1st.2nd$conc.1st.540))
setwd("/user/xinyugao/ELISA_project/")
pdf("compare_concentration_450-540.pdf")
plot(res.1st.2nd$conc.1st.540, res.1st.2nd$conc.2nd.540, pch=19, col="black",
     xlab="1st Concentration", ylab="2nd Concentration",
     main = "scatterplot of 2nd and 1st concentrations (450-540nm)")
abline(lm(res.1st.2nd$conc.2nd.540 ~ res.1st.2nd$conc.1st.540), col = "red")
dev.off()

# scatterplot of 570
summary(lm(res.1st.2nd$conc.2nd.570 ~ res.1st.2nd$conc.1st.570))
setwd("/user/xinyugao/ELISA_project/")
pdf("compare_concentration_450-570.pdf")
plot(res.1st.2nd$conc.1st.570, res.1st.2nd$conc.2nd.570, pch=19, col="black",
     xlab="1st Concentration", ylab="2nd Concentration",
     main = "scatterplot of 2nd and 1st concentrations (450-570nm)")
abline(lm(res.1st.2nd$conc.2nd.570 ~ res.1st.2nd$conc.1st.570), col = "red")
dev.off()







