rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)
# library(hrbrthemes)

filepath <- "/projects/rpci/qzhu/users-qzhu/grant/UACA_R01/ELISA/"
PT.genotype <- read.table(paste0(filepath, "5th_ELISA.select_samples.txt"), header = TRUE, sep = "\t")
PT.genotype.rm <- PT.genotype[- which(PT.genotype$CollectionID=="Cl-00044311" & PT.genotype$COLLECTIONBARCODE == "1002201"), ]

out.PT.genotype <- PT.genotype.rm %>% select(SubjectID, rs6494889)
dat.450_570 <- read.table(paste0(filepath, "5thELISA_450nm-570nm"), header = TRUE, colClasses = c("numeric", "numeric", "numeric"))
dat.450_540 <- read.table(paste0(filepath, "5thELISA_450nm-540nm"), header = TRUE, colClasses = c("numeric", "numeric", "numeric"))


ID.info <- readxl::read_excel("/projects/rpci/qzhu/users-qzhu/grant/UACA_R01/ELISA/DistSamplesList_nowaka_Sep-16-2022-11-25-05.xls") # collection_alias
ID.info <- ID.info %>% select(Subjectid, `Collection Alias`, Container) 
colnames(ID.info) <- c("Subjectid", "ID_alias", "Container")
ID.info <- unique( ID.info )
ID.info$ID_alias <- as.numeric(ID.info$ID_alias)

dat.450_540 <- left_join(dat.450_540, ID.info, by=c("ID"="ID_alias"))  
dat.450_570 <- left_join(dat.450_570, ID.info, by=c("ID"="ID_alias"))
length(unique(dat.450_570$Subjectid)) # 39 PTs
length(unique(dat.450_570$ID)) # 40 records


dat.450_540.genotype <- inner_join(dat.450_540, out.PT.genotype, by=c("Subjectid"="SubjectID"))
dat.450_570.genotype <- inner_join(dat.450_570, out.PT.genotype, by=c("Subjectid"="SubjectID"))
length(unique(dat.450_540.genotype$Subjectid))  # 38
length(unique(dat.450_570.genotype$Subjectid))  # 38 

dat.450_540.genotype$cleaved_par4 <- dat.450_540.genotype$original - 0.5*dat.450_540.genotype$X30K
dat.450_570.genotype$cleaved_par4 <- dat.450_570.genotype$original - 0.5*dat.450_570.genotype$X30K


# Please test the original, 30K, and cleaved Par4 concentration respectively for normal distribution
# Shapiro-Wilk's test: The null hypothesis of these tests is that sample distribution is normal. If the test is significant, the distribution is non-normal.

# 450-540 original
shapiro.test(dat.450_540.genotype$original)$p.value
# p-value= 0.005883792 < 0.05, we reject the normality assumption
hist(dat.450_540.genotype$original, br=10, xlab="concentration", main="original (450-540nm), normality test p-value=0.0059")

# 450-540 30K
shapiro.test(dat.450_540.genotype$X30K)$p.value 
# p-value=0.01767933 < 0.05, we reject the normality assumption
hist(dat.450_540.genotype$X30K, br=10, xlab="concentration", main="30K (450-540nm), normality test p-value=0.0177")

# 450-540 cleaved
shapiro.test(dat.450_540.genotype$cleaved_par4)$p.value 
# p-value=0.5443517 > 0.05, do not reject H0, we can assume the normality
hist(dat.450_540.genotype$cleaved_par4, br=10, xlab="concentration", main="cleaved_par4 (450-540nm), normality test p-value=0.5443")




dat.hist.450_540 <- dat.450_540.genotype %>% select(original, X30K, cleaved_par4) %>%
  gather(key="text", value="value")
# plot
setwd("/user/xinyugao/ELISA_project/")
pdf(paste0("5th_hist_450_540nm.pdf"))
dat.hist.450_540 %>%
  ggplot( aes(x=value, fill=text)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#69b3a2", "#404080", "goldenrod2")) +
  xlab("Concentrations (450-540nm)") +
  ylab("Frequency") +
  labs(fill="")
dev.off()

dat.hist.450_570 <- dat.450_570.genotype %>% select(original, X30K, cleaved_par4) %>%
  gather(key="text", value="value")
# plot
setwd("/user/xinyugao/ELISA_project/")
pdf(paste0("5th_hist_450_570nm.pdf"))
dat.hist.450_570 %>%
  ggplot( aes(x=value, fill=text)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#69b3a2", "#404080", "goldenrod2")) +
  xlab("Concentrations (450-570nm)") +
  ylab("Frequency") +
  labs(fill="")
dev.off()


# 450-570 original
shapiro.test(dat.450_570.genotype$original)$p.value
# p-value= 0.005724294 < 0.05, we reject the normality assumption
hist(dat.450_570.genotype$original, br=10, xlab="concentration", main="original (450-570nm), normality test p-value=0.0057")

# 450-570 30K
shapiro.test(dat.450_570.genotype$X30K)$p.value 
# p-value=0.01729104 < 0.05, we reject the normality assumption
hist(dat.450_570.genotype$X30K, br=10, xlab="concentration", main="30K (450-570nm), normality test p-value=0.0173")

# 450-570 cleaved
shapiro.test(dat.450_570.genotype$cleaved_par4)$p.value 
# p-value=0.5425903 > 0.05, do not reject H0, we can assume the normality
hist(dat.450_570.genotype$cleaved_par4, br=10, xlab="concentration", main="cleaved_par4 (450-570nm), normality test p-value=0.543")

## log-transform?

# 450-540nm
shapiro.test(log(dat.450_540.genotype$original))$p.value
# p-value=0.07378116 > 0.05, do not reject H0, we can assume the normality
hist(log(dat.450_540.genotype$original), br=10, xlab="log(concentration)", main="log(original) (450-570nm), normality test p-value=0.074")

shapiro.test(log(dat.450_540.genotype$X30K))$p.value
# p-value=0.1427699 > 0.05, do not reject H0, we can assume the normality

# 450-570nm
shapiro.test(log(dat.450_570.genotype$original))$p.value
# p-value=0.1038735 > 0.05, do not reject H0, we can assume the normality
shapiro.test(log(dat.450_570.genotype$X30K))$p.value
# p-value=0.1925116 > 0.05, do not reject H0, we can assume the normality


# please also generate boxplot of each of the 3 concentrations by the 3 genotype groups with 
# t-test and wilcox test pvalues for comparision between genotype 1 and 0 as well as between genotype 2 and 0.

# table(dat.450_540.genotype$rs6494889)
# 0  1  2 
# 19 17  2 
# table(dat.450_570.genotype$rs6494889)
# 0  1  2 
# 19 17  2

###################################################################################################


get_test_res <- function(data, container) {
  #---------------------- 1 vs 0 ---------------------
  x = data$original[(data$rs6494889==0) & (data$Container == container)]
  y = data$original[(data$rs6494889==1) & (data$Container == container)]
  t.test_0vs1.original  = t.test(x, y)$p.value
  wilcox_0vs1.original  = wilcox.test(x, y)$p.value
  
  x = data$X30K[(data$rs6494889==0) & (data$Container == container)]
  y = data$X30K[(data$rs6494889==1) & (data$Container == container)]
  t.test_0vs1.30K  = t.test(x, y)$p.value
  wilcox_0vs1.30K  = wilcox.test(x, y)$p.value
  
  x = data$cleaved_par4[(data$rs6494889==0) & (data$Container == container)]
  y = data$cleaved_par4[(data$rs6494889==1) & (data$Container == container)]
  t.test_0vs1.cleaved_par4  = t.test(x, y)$p.value
  wilcox_0vs1.cleaved_par4  = wilcox.test(x, y)$p.value
  #---------------------- 1 vs 0 ---------------------
  
  #---------------------- 1+2 vs 0 ---------------------
  x = data$original[(data$rs6494889==0) & (data$Container == container)]
  y = data$original[(data$rs6494889!=0) & (data$Container == container)]
  t.test_0vs12.original  = t.test(x, y)$p.value
  wilcox_0vs12.original  = wilcox.test(x, y)$p.value
  
  x = data$X30K[(data$rs6494889==0) & (data$Container == container)]
  y = data$X30K[(data$rs6494889!=0) & (data$Container == container)]
  t.test_0vs12.30K  = t.test(x, y)$p.value
  wilcox_0vs12.30K  = wilcox.test(x, y)$p.value
  
  x = data$cleaved_par4[(data$rs6494889==0) & (data$Container == container)]
  y = data$cleaved_par4[(data$rs6494889!=0) & (data$Container == container)]
  t.test_0vs12.cleaved_par4  = t.test(x, y)$p.value
  wilcox_0vs12.cleaved_par4  = wilcox.test(x, y)$p.value
  #---------------------- 1+2 vs 0 ---------------------
  options(digits=2)
  return(list("t.test_0vs1.original"=t.test_0vs1.original, "wilcox_0vs1.original"=wilcox_0vs1.original,
              "t.test_0vs1.30K "=t.test_0vs1.30K , "wilcox_0vs1.30K"=wilcox_0vs1.30K,
              "t.test_0vs1.cleaved_par4"=t.test_0vs1.cleaved_par4, "wilcox_0vs1.cleaved_par4"=wilcox_0vs1.cleaved_par4,
              "t.test_0vs12.original"=t.test_0vs12.original, "wilcox_0vs12.original"=wilcox_0vs12.original,
              "t.test_0vs12.30K "=t.test_0vs12.30K , "wilcox_0vs12.30K"=wilcox_0vs12.30K,
              "t.test_0vs12.cleaved_par4"=t.test_0vs12.cleaved_par4, "wilcox_0vs12.cleaved_par4"=wilcox_0vs12.cleaved_par4))
}


res.green <- get_test_res(data = dat.450_540.genotype, container = "CBS Straw Green")
res.white <- get_test_res(data = dat.450_540.genotype, container = "CBS Straw White")

dat.450_540.genotype$rs6494889 <- as.factor(dat.450_540.genotype$rs6494889)
#----------------------------------------------------------------------------------------------
# original
anno <- data.frame(xstar = c(2, 2), ystar = c(175, 175),
                   lab = c(paste0("'1' vs '0', t-test: p=",round(res.green$t.test_0vs1.original,2) ,
                                  "; wilcox: p=", round(res.green$wilcox_0vs1.original, 2), "\n", 
                                  "'1+2' vs '0', t-test: p=",round(res.green$t.test_0vs12.original,2) ,
                                  "; wilcox: p=", round(res.green$wilcox_0vs12.original,2)), 
                           paste0("'1' vs '0', t-test: p=",round(res.white$t.test_0vs1.original,2) ,
                                  "; wilcox: p=", round(res.white$wilcox_0vs1.original, 2), "\n", 
                                  "'1+2' vs '0', t-test: p=",round(res.white$t.test_0vs12.original,2) ,
                                  "; wilcox: p=", round(res.white$wilcox_0vs12.original,2))
                           ),
                   Container = names(table(dat.450_540.genotype$Container)))

setwd("/user/xinyugao/ELISA_project/")
pdf(paste0("5th.rm_Sample_boxplot_ContainersANDGenotypes_450_540.original.pdf"), height=4)
ggplot(dat.450_540.genotype, aes(x = rs6494889, y = original)) +
  geom_boxplot(aes(fill=rs6494889), alpha=0.3, outlier.shape = NA) + geom_jitter(aes(color=rs6494889), show.legend = F)+ 
  xlab("Genotypes") + ylab("Concentration") +  labs(fill='Genotypes')  + 
  scale_y_continuous(limits=c(0,175)) +
  theme(legend.position="right") +theme(legend.title = element_text(size=10,  face="bold")) +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        strip.background = element_rect(fill = "white", colour = "black"),
        legend.key = element_blank()) + 
  geom_text(data = anno, aes(x = xstar,  y = ystar, label = lab), size=2.5) +
  facet_wrap(.~Container, nrow=1)  +   theme(panel.spacing = unit(1, "lines")) 
dev.off()
#---------------------------------------------------------------------------------------------

# 30K
anno <- data.frame(xstar = c(2, 2), ystar = c(175, 175),
                   lab = c(paste0("'1' vs '0', t-test: p=",round(res.green$t.test_0vs1.30K ,2) ,
                                  "; wilcox: p=", round(res.green$wilcox_0vs1.30K, 2), "\n", 
                                  "'1+2' vs '0', t-test: p=",round(res.green$t.test_0vs12.30K,2) ,
                                  "; wilcox: p=", round(res.green$wilcox_0vs12.30K,2)), 
                           paste0("'1' vs '0', t-test: p=",round(res.white$t.test_0vs1.30K,2) ,
                                  "; wilcox: p=", round(res.white$wilcox_0vs1.30K, 2), "\n", 
                                  "'1+2' vs '0', t-test: p=",round(res.white$t.test_0vs12.30K,2) ,
                                  "; wilcox: p=", round(res.white$wilcox_0vs12.30K,2))
                   ),
                   Container = names(table(dat.450_540.genotype$Container)))

setwd("/user/xinyugao/ELISA_project/")
pdf(paste0("5th.rm_Sample_boxplot_ContainersANDGenotypes_450_540.30K.pdf"), height=4)
ggplot(dat.450_540.genotype, aes(x = rs6494889, y = X30K)) +
  geom_boxplot(aes(fill=rs6494889), alpha=0.3, outlier.shape = NA) + geom_jitter(aes(color=rs6494889), show.legend = F)+ 
  xlab("Genotypes") + ylab("Concentration") +  labs(fill='Genotypes')  + 
  scale_y_continuous(limits=c(0,205)) +
  theme(legend.position="right") +theme(legend.title = element_text(size=10,  face="bold")) +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        strip.background = element_rect(fill = "white", colour = "black"),
        legend.key = element_blank()) + 
  geom_text(data = anno, aes(x = xstar,  y = ystar, label = lab), size=2.5) +
  facet_wrap(.~Container, nrow=1)  +   theme(panel.spacing = unit(1, "lines")) 
dev.off()
#---------------------------------------------------------------------------------------------

# cleaved_par4
anno <- data.frame(xstar = c(2, 2), ystar = c(135,135),
                   lab = c(paste0("'1' vs '0', t-test: p=",round(res.green$t.test_0vs1.cleaved_par4 ,2) ,
                                  "; wilcox: p=", round(res.green$wilcox_0vs1.cleaved_par4, 2), "\n", 
                                  "'1+2' vs '0', t-test: p=",round(res.green$t.test_0vs12.cleaved_par4,2) ,
                                  "; wilcox: p=", round(res.green$wilcox_0vs12.cleaved_par4,2)), 
                           paste0("'1' vs '0', t-test: p=",round(res.white$t.test_0vs1.cleaved_par4,2) ,
                                  "; wilcox: p=", round(res.white$wilcox_0vs1.cleaved_par4, 2), "\n", 
                                  "'1+2' vs '0', t-test: p=",round(res.white$t.test_0vs12.cleaved_par4,2) ,
                                  "; wilcox: p=", round(res.white$wilcox_0vs12.cleaved_par4,2))
                   ),
                   Container = names(table(dat.450_540.genotype$Container)))

setwd("/user/xinyugao/ELISA_project/")
pdf(paste0("5th.rm_Sample_boxplot_ContainersANDGenotypes_450_540.cleaved_par4.pdf"), height=4)
ggplot(dat.450_540.genotype, aes(x = rs6494889, y = cleaved_par4)) +
  geom_boxplot(aes(fill=rs6494889), alpha=0.3, outlier.shape = NA) + 
  geom_jitter(aes(color=rs6494889), show.legend = F)+ 
  xlab("Genotypes") + ylab("Concentration") +  labs(fill='Genotypes')  + 
  scale_y_continuous(limits=c(-85,135)) +
  theme(legend.position="right") +theme(legend.title = element_text(size=10,  face="bold")) +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        strip.background = element_rect(fill = "white", colour = "black"),
        legend.key = element_blank()) + 
  geom_text(data = anno, aes(x = xstar,  y = ystar, label = lab), size=2.5) +
  facet_wrap(.~Container, nrow=1)  +   theme(panel.spacing = unit(1, "lines")) 
dev.off()


###################################################################################################





gen.dat.boxplot <-function(data) {
  names <- c(rep("genotype=0", 19) , rep("genotype=1", 17) , rep("genotype=2", 2))
  value.original <- c( data$original[data$rs6494889==0] , 
                       data$original[data$rs6494889==1],
                       data$original[data$rs6494889==2])
  value.30K <- c( data$X30K[data$rs6494889==0] , 
                  data$X30K[data$rs6494889==1],
                  data$X30K[data$rs6494889==2])
  value.cleaved <- c( data$cleaved_par4[data$rs6494889==0] , 
                      data$cleaved_par4[data$rs6494889==1],
                      data$cleaved_par4[data$rs6494889==2])
  df <- data.frame("names"=names, "original"=value.original, "30K"=value.30K, "cleaved"=value.cleaved)
  
}

dat.boxplot.450_540 <- gen.dat.boxplot(dat.450_540.genotype)



setwd("/user/xinyugao/ELISA_project/")
pdf(paste0("5th.rm_Sample_boxplot_genotypes_450_540.original.pdf"))
data <- dat.boxplot.450_540
boxplot(data$original ~ data$names , col=terrain.colors(4) , ylim=c(0, max(data$original) + 50),outline=FALSE,
        xlab = "Genotype Groups", ylab = "Concentration")
# Add data points
mylevels <- levels(factor(data$names))
levelProportions <- table(data$names)/nrow(data)
for(i in 1:length(mylevels)){
  
  thislevel <- mylevels[i]
  thisvalues <- data$original[data$names==thislevel]
  
  # take the x-axis indices and add a jitter, proportional to the N in each level
  myjitter <- jitter(rep(i, length(thisvalues)), amount=levelProportions[i]/2)
  points(myjitter, thisvalues, pch=20, col=rgb(0,0,0,.9)) 
}
legend("topright",pt.cex = 2, border = F,
       text.col = "black", 
       horiz = F , 
       legend = c("genotype=0 (n=19)", "genotype=1 (n=17)", "genotype=2 (n=2)"),
       fill = terrain.colors(4))
text(1, 190, paste0("'1' vs '0', t-test: p=", round(t.test(data$original[data$names=="genotype=0"], data$original[data$names=="genotype=1"])$p.value, 2),
                    "; wilcox: p=", round(wilcox.test(data$original[data$names=="genotype=0"], data$original[data$names=="genotype=1"])$p.value, 2)),  cex = 0.7)
text(1, 200, paste0("'2' vs '0', t-test: p=", round(t.test(data$original[data$names=="genotype=0"], data$original[data$names=="genotype=2"])$p.value, 2),
                    "; wilcox: p=", round(wilcox.test(data$original[data$names=="genotype=0"], data$original[data$names=="genotype=2"])$p.value,2)),  cex = 0.7)
text(1, 210, paste0("'1+2' vs '0', t-test: p=", round(t.test(data$original[data$names=="genotype=0"], data$original[data$names != "genotype=0"])$p.value, 2),
                    "; wilcox: p=", round(wilcox.test(data$original[data$names=="genotype=0"], data$original[data$names != "genotype=0"])$p.value, 2)),  cex = 0.7)

dev.off()


setwd("/user/xinyugao/ELISA_project/")
pdf(paste0("5th.rm_Sample_boxplot_genotypes_450_540.30K.pdf"))
data <- dat.boxplot.450_540
boxplot(data$X30K ~ data$names , col=terrain.colors(4) , ylim=c(0, max(data$X30K) + 20),outline=FALSE,
        xlab = "Genotype Groups", ylab = "Concentration")
# Add data points
mylevels <- levels(factor(data$names))
levelProportions <- table(data$names)/nrow(data)
for(i in 1:length(mylevels)){
  
  thislevel <- mylevels[i]
  thisvalues <- data$X30K[data$names==thislevel]
  
  # take the x-axis indices and add a jitter, proportional to the N in each level
  myjitter <- jitter(rep(i, length(thisvalues)), amount=levelProportions[i]/2)
  points(myjitter, thisvalues, pch=20, col=rgb(0,0,0,.9)) 
}
legend("topright",pt.cex = 2, border = F,
       text.col = "black", 
       horiz = F , 
       legend = c("genotype=0 (n=19)", "genotype=1 (n=17)", "genotype=2 (n=2)"),
       fill = terrain.colors(4))
text(1, 190, paste0("'1' vs '0', t-test: p=", round(t.test(data$X30K[data$names=="genotype=0"], data$X30K[data$names=="genotype=1"])$p.value, 2),
                    "; wilcox: p=", round(wilcox.test(data$X30K[data$names=="genotype=0"], data$X30K[data$names=="genotype=1"])$p.value, 2)), cex = 0.7)
text(1, 200, paste0("'2' vs '0', t-test: p=", round(t.test(data$X30K[data$names=="genotype=0"], data$X30K[data$names=="genotype=2"])$p.value, 2),
                    "; wilcox: p=", round(wilcox.test(data$X30K[data$names=="genotype=0"], data$X30K[data$names=="genotype=2"])$p.value, 2)), cex = 0.7)
text(1, 210, paste0("'1+2' vs '0', t-test: p=", round(t.test(data$X30K[data$names=="genotype=0"], data$X30K[data$names!="genotype=0"])$p.value,2),
                    "; wilcox: p=", round(wilcox.test(data$X30K[data$names=="genotype=0"], data$X30K[data$names!="genotype=0"])$p.value, 2)), cex = 0.7)

dev.off()



setwd("/user/xinyugao/ELISA_project/")
pdf(paste0("5th.rm_Sample_boxplot_genotypes_450_540.cleaved_par4.pdf"))
data <- dat.boxplot.450_540
boxplot(data$cleaved~ data$names , col=terrain.colors(4) , ylim=c(min(data$cleaved), max(data$cleaved) + 50),outline=FALSE,
        xlab = "Genotype Groups", ylab = "Concentration")
# Add data points
mylevels <- levels(factor(data$names))
levelProportions <- table(data$names)/nrow(data)
for(i in 1:length(mylevels)){
  
  thislevel <- mylevels[i]
  thisvalues <- data$cleaved[data$names==thislevel]
  
  # take the x-axis indices and add a jitter, proportional to the N in each level
  myjitter <- jitter(rep(i, length(thisvalues)), amount=levelProportions[i]/2)
  points(myjitter, thisvalues, pch=20, col=rgb(0,0,0,.9)) 
}
legend("topright",pt.cex = 2, border = F,
       text.col = "black", 
       horiz = F , 
       legend = c("genotype=0 (n=19)", "genotype=1 (n=17)", "genotype=2 (n=2)"),
       fill = terrain.colors(4))
text(1, 140, paste0("'1' vs '0', t-test: p=", round(t.test(data$cleaved[data$names=="genotype=0"], data$cleaved[data$names=="genotype=1"])$p.value, 2),
                    "; wilcox: p=", round(wilcox.test(data$cleaved[data$names=="genotype=0"], data$cleaved[data$names=="genotype=1"])$p.value, 2)), cex = 0.7)
text(1, 150, paste0("'2' vs '0', t-test: p=", round(t.test(data$cleaved[data$names=="genotype=0"], data$cleaved[data$names=="genotype=2"])$p.value, 2),
                    "; wilcox: p=", round(wilcox.test(data$cleaved[data$names=="genotype=0"], data$cleaved[data$names=="genotype=2"])$p.value, 2)), cex = 0.7)
text(1, 160, paste0("'1+2' vs '0', t-test: p=", round(t.test(data$cleaved[data$names=="genotype=0"], data$cleaved[data$names!="genotype=0"])$p.value, 2),
                    "; wilcox: p=", round(wilcox.test(data$cleaved[data$names=="genotype=0"], data$cleaved[data$names!="genotype=0"])$p.value, 2)), cex = 0.7)

dev.off()


#--------------------------------------------------------------------------------------
dat.boxplot.450_570 <- gen.dat.boxplot(dat.450_570.genotype)

setwd("/user/xinyugao/ELISA_project/")
pdf(paste0("5th.rm_Sample_boxplot_genotypes_450_570.original.pdf"))
data <- dat.boxplot.450_570
boxplot(data$original ~ data$names , col=terrain.colors(4) , ylim=c(0, max(data$original) + 50),outline=FALSE,
        xlab = "Genotype Groups", ylab = "Concentration")
# Add data points
mylevels <- levels(factor(data$names))
levelProportions <- table(data$names)/nrow(data)
for(i in 1:length(mylevels)){
  
  thislevel <- mylevels[i]
  thisvalues <- data$original[data$names==thislevel]
  
  # take the x-axis indices and add a jitter, proportional to the N in each level
  myjitter <- jitter(rep(i, length(thisvalues)), amount=levelProportions[i]/2)
  points(myjitter, thisvalues, pch=20, col=rgb(0,0,0,.9)) 
}
legend("topright",pt.cex = 2, border = F,
       text.col = "black", 
       horiz = F , 
       legend = c("genotype=0 (n=19)", "genotype=1 (n=17)", "genotype=2 (n=2)"),
       fill = terrain.colors(4))
text(1, 190, paste0("'1' vs '0', t-test: p=", round(t.test(data$original[data$names=="genotype=0"], data$original[data$names=="genotype=1"])$p.value, 2),
                    "; wilcox: p=", round(wilcox.test(data$original[data$names=="genotype=0"], data$original[data$names=="genotype=1"])$p.value, 2)),  cex = 0.7)
text(1, 200, paste0("'2' vs '0', t-test: p=", round(t.test(data$original[data$names=="genotype=0"], data$original[data$names=="genotype=2"])$p.value, 2),
                    "; wilcox: p=", round(wilcox.test(data$original[data$names=="genotype=0"], data$original[data$names=="genotype=2"])$p.value, 2)),  cex = 0.7)
text(1, 210, paste0("'1+2' vs '0', t-test: p=", round(t.test(data$original[data$names=="genotype=0"], data$original[data$names!="genotype=0"])$p.value, 2),
                    "; wilcox: p=", round(wilcox.test(data$original[data$names=="genotype=0"], data$original[data$names!="genotype=0"])$p.value, 2)),  cex = 0.7)
dev.off()



setwd("/user/xinyugao/ELISA_project/")
pdf(paste0("5th.rm_Sample_boxplot_genotypes_450_570.30K.pdf"))
data <- dat.boxplot.450_570
boxplot(data$X30K ~ data$names , col=terrain.colors(4) , ylim=c(0, max(data$X30K) + 20),outline=FALSE,
        xlab = "Genotype Groups", ylab = "Concentration")
# Add data points
mylevels <- levels(factor(data$names))
levelProportions <- table(data$names)/nrow(data)
for(i in 1:length(mylevels)){
  
  thislevel <- mylevels[i]
  thisvalues <- data$X30K[data$names==thislevel]
  
  # take the x-axis indices and add a jitter, proportional to the N in each level
  myjitter <- jitter(rep(i, length(thisvalues)), amount=levelProportions[i]/2)
  points(myjitter, thisvalues, pch=20, col=rgb(0,0,0,.9)) 
}
legend("topright",pt.cex = 2, border = F,
       text.col = "black", 
       horiz = F , 
       legend = c("genotype=0 (n=19)", "genotype=1 (n=17)", "genotype=2 (n=2)"),
       fill = terrain.colors(4))
text(1, 190, paste0("'1' vs '0', t-test: p=", round(t.test(data$X30K[data$names=="genotype=0"], data$X30K[data$names=="genotype=1"])$p.value, 2),
                    "; wilcox: p=", round(wilcox.test(data$X30K[data$names=="genotype=0"], data$X30K[data$names=="genotype=1"])$p.value, 2)), cex = 0.7)
text(1, 200, paste0("'2' vs '0', t-test: p=", round(t.test(data$X30K[data$names=="genotype=0"], data$X30K[data$names=="genotype=2"])$p.value, 2),
                    "; wilcox: p=", round(wilcox.test(data$X30K[data$names=="genotype=0"], data$X30K[data$names=="genotype=2"])$p.value, 2)), cex = 0.7)
text(1, 210, paste0("'1+2' vs '0', t-test: p=", round(t.test(data$X30K[data$names=="genotype=0"], data$X30K[data$names!="genotype=0"])$p.value, 2),
                    "; wilcox: p=", round(wilcox.test(data$X30K[data$names=="genotype=0"], data$X30K[data$names!="genotype=0"])$p.value, 2)), cex = 0.7)
dev.off()


setwd("/user/xinyugao/ELISA_project/")
pdf(paste0("5th.rm_Sample_boxplot_genotypes_450_570.cleaved_par4.pdf"))
data <- dat.boxplot.450_570
boxplot(data$cleaved~ data$names , col=terrain.colors(4) , ylim=c(min(data$cleaved), max(data$cleaved) + 50),outline=FALSE,
        xlab = "Genotype Groups", ylab = "Concentration")
# Add data points
mylevels <- levels(factor(data$names))
levelProportions <- table(data$names)/nrow(data)
for(i in 1:length(mylevels)){
  
  thislevel <- mylevels[i]
  thisvalues <- data$cleaved[data$names==thislevel]
  
  # take the x-axis indices and add a jitter, proportional to the N in each level
  myjitter <- jitter(rep(i, length(thisvalues)), amount=levelProportions[i]/2)
  points(myjitter, thisvalues, pch=20, col=rgb(0,0,0,.9)) 
}
legend("topright",pt.cex = 2, border = F,
       text.col = "black", 
       horiz = F , 
       legend = c("genotype=0 (n=19)", "genotype=1 (n=17)", "genotype=2 (n=2)"),
       fill = terrain.colors(4))
text(1, 140, paste0("'1' vs '0', t-test: p=", round(t.test(data$cleaved[data$names=="genotype=0"], data$cleaved[data$names=="genotype=1"])$p.value, 2),
                    "; wilcox: p=", round(wilcox.test(data$cleaved[data$names=="genotype=0"], data$cleaved[data$names=="genotype=1"])$p.value, 2)), cex = 0.7)
text(1, 150, paste0("'2' vs '0', t-test: p=", round(t.test(data$cleaved[data$names=="genotype=0"], data$cleaved[data$names=="genotype=2"])$p.value, 2),
                    "; wilcox: p=", round(wilcox.test(data$cleaved[data$names=="genotype=0"], data$cleaved[data$names=="genotype=2"])$p.value, 2)), cex = 0.7)
text(1, 160, paste0("'1+2' vs '0', t-test: p=", round(t.test(data$cleaved[data$names=="genotype=0"], data$cleaved[data$names!="genotype=0"])$p.value, 2),
                    "; wilcox: p=", round(wilcox.test(data$cleaved[data$names=="genotype=0"], data$cleaved[data$names!="genotype=0"])$p.value, 2)), cex = 0.7)

dev.off()




# 1. please run linear regression between Y=original concentration and X=0.5*(30K concentration) 
# as well as pearson and spearman correlation coefficients for all samples, genotype=0, genotype=1 samples respectively after removing that one sample.
# 2. please linear regression between Y=concentration of cleaved Par4 and X=0.5*(30K concentration) 
# as well as pearson and spearman correlation coefficients for all samples, genotype=0, genotype=1 samples respectively after removing that one sample.




get_LRs <- function(input, container) {
  #----------------------------------------#
  # Y=original concentration and X=0.5*(30K concentration) 
  coef.pearson.original_30K <- c(); coef.spearman.original_30K <- c()
  ## all samples
  y <- input$original[input$Container == container]; x <- input$X30K[input$Container == container] / 2
  LR.sum <- summary(lm(y ~ x))
  LR.out.all <- LR.sum$coefficients
  coef.pearson.all <- cor(x, y, method = c("pearson"))
  coef.spearman.all <- cor(x, y, method = c("spearman"))

  ## genotype=0
  y <- input$original[input$rs6494889==0 & input$Container == container]; 
  x <- input$X30K[input$rs6494889==0 & input$Container == container] / 2
  LR.sum <- summary(lm(y ~ x))
  LR.out.0 <- LR.sum$coefficients
  coef.pearson.0 <- cor(x, y, method = c("pearson"))
  coef.spearman.0 <- cor(x, y, method = c("spearman"))
  
  ## genotype=1
  y <- input$original[input$rs6494889==1 & input$Container == container]; 
  x <- input$X30K[input$rs6494889==1 & input$Container == container] / 2
  LR.sum <- summary(lm(y ~ x))
  LR.out.1 <- LR.sum$coefficients
  coef.pearson.1 <- cor(x, y, method = c("pearson"))
  coef.spearman.1 <- cor(x, y, method = c("spearman"))
  
  ## genotype=1+2
  y <- input$original[input$rs6494889!=0 & input$Container == container]; 
  x <- input$X30K[input$rs6494889!=0 & input$Container == container] / 2
  LR.sum <- summary(lm(y ~ x))
  LR.out.12 <- LR.sum$coefficients
  coef.pearson.12 <- cor(x, y, method = c("pearson"))
  coef.spearman.12 <- cor(x, y, method = c("spearman"))
  
  
  # output
  LR.original_30K <- rbind(LR.out.all, LR.out.0, LR.out.1, LR.out.12)
  coef.pearson.original_30K <- c(coef.pearson.all, coef.pearson.0, coef.pearson.1, coef.pearson.12)
  coef.spearman.original_30K <- c(coef.spearman.all,  coef.spearman.0, coef.spearman.1, coef.spearman.12)
  
  # Y=concentration of cleaved Par4 and X=0.5*(30K concentration) 
  coef.pearson.cleaved_30K <- c(); coef.spearman.cleaved_30K <- c()
  ## all samples
  y <- input$cleaved_par4[input$Container == container]; x <- input$X30K[input$Container == container] / 2
  LR.sum <- summary(lm(y ~ x))
  LR.out.all <- LR.sum$coefficients
  coef.pearson.all <- cor(x, y, method = c("pearson"))
  coef.spearman.all <- cor(x, y, method = c("spearman"))

  ## genotype=0
  y <- input$cleaved_par4[input$rs6494889==0 & input$Container == container]; 
  x <- input$X30K[input$rs6494889==0 & input$Container == container] / 2
  LR.sum <- summary(lm(y ~ x))
  LR.out.0 <- LR.sum$coefficients
  coef.pearson.0 <- cor(x, y, method = c("pearson"))
  coef.spearman.0 <- cor(x, y, method = c("spearman"))
  
  # plot(x, y, xlab = "0.5*30K", ylab = "cleaved_par4")
  # abline(lm(y ~ x), col = "red")
  
  
  ## genotype=1
  y <- input$cleaved_par4[input$rs6494889==1 & input$Container == container]; 
  x <- input$X30K[input$rs6494889==1 & input$Container == container] / 2
  LR.sum <- summary(lm(y ~ x))
  LR.out.1 <- LR.sum$coefficients
  coef.pearson.1 <- cor(x, y, method = c("pearson"))
  coef.spearman.1 <- cor(x, y, method = c("spearman"))
  
  ## genotype=1+2
  y <- input$cleaved_par4[input$rs6494889!=0 & input$Container == container]; 
  x <- input$X30K[input$rs6494889!=0 & input$Container == container] / 2
  LR.sum <- summary(lm(y ~ x))
  LR.out.12 <- LR.sum$coefficients
  coef.pearson.12 <- cor(x, y, method = c("pearson"))
  coef.spearman.12 <- cor(x, y, method = c("spearman"))
  
  # output
  LR.cleaved_30K <- rbind(LR.out.all, LR.out.0, LR.out.1, LR.out.12)
  coef.pearson.cleaved_30K <- c(coef.pearson.all, coef.pearson.0, coef.pearson.1, coef.pearson.12)
  coef.spearman.cleaved_30K <- c(coef.spearman.all,  coef.spearman.0, coef.spearman.1, coef.spearman.12)
  
  return(list("LR.original_30K"=LR.original_30K, 
              "LR.cleaved_30K"=LR.cleaved_30K, 
              "ceof"=data.frame(coef.pearson.original_30K, coef.spearman.original_30K,
                                coef.pearson.cleaved_30K, coef.spearman.cleaved_30K)))
  
}


options(digits = 2)
get_LRs(input=dat.450_540.genotype, container = "CBS Straw Green")
get_LRs(input=dat.450_540.genotype, container = "CBS Straw White")

get_LRs(input=dat.450_570.genotype, container = "CBS Straw Green")
get_LRs(input=dat.450_570.genotype, container = "CBS Straw White")


# could you please also test linear regression between Y=original concentration and X=cleaved par=4?


get_LRs <- function(input, container) {
  #----------------------------------------#
  # Y=original concentration and X=cleaved par=4
  coef.pearson <- c(); coef.spearman <- c()
  ## all samples
  y <- input$original[input$Container == container]; x <- input$cleaved_par4[input$Container == container]
  LR.sum <- summary(lm(y ~ x))
  LR.out.all <- LR.sum$coefficients
  coef.pearson.all <- cor(x, y, method = c("pearson"))
  coef.spearman.all <- cor(x, y, method = c("spearman"))
  
  ## genotype=0
  y <- input$original[input$rs6494889==0 & input$Container == container]; 
  x <- input$cleaved_par4[input$rs6494889==0 & input$Container == container] 
  LR.sum <- summary(lm(y ~ x))
  LR.out.0 <- LR.sum$coefficients
  coef.pearson.0 <- cor(x, y, method = c("pearson"))
  coef.spearman.0 <- cor(x, y, method = c("spearman"))
  
  ## genotype=1
  y <- input$original[input$rs6494889==1 & input$Container == container]; 
  x <- input$cleaved_par4[input$rs6494889==1 & input$Container == container] 
  LR.sum <- summary(lm(y ~ x))
  LR.out.1 <- LR.sum$coefficients
  coef.pearson.1 <- cor(x, y, method = c("pearson"))
  coef.spearman.1 <- cor(x, y, method = c("spearman"))
  
  ## genotype=1+2
  y <- input$original[input$rs6494889!=0 & input$Container == container]; 
  x <- input$cleaved_par4[input$rs6494889!=0 & input$Container == container] 
  LR.sum <- summary(lm(y ~ x))
  LR.out.12 <- LR.sum$coefficients
  coef.pearson.12 <- cor(x, y, method = c("pearson"))
  coef.spearman.12 <- cor(x, y, method = c("spearman"))
  
  # output
  LR <- rbind(LR.out.all, LR.out.0, LR.out.1, LR.out.12)
  coef.pearson <- c(coef.pearson.all, coef.pearson.0, coef.pearson.1, coef.pearson.12)
  coef.spearman <- c(coef.spearman.all,  coef.spearman.0, coef.spearman.1, coef.spearman.12)
  
  return(list("LR"=LR, 
              "coef.pearson"=data.frame("name"=c("all", "0", "1", "1+2"), "coef"=coef.pearson),
              "coef.spearman"=data.frame("name"=c("all", "0", "1", "1+2"), "coef"=coef.spearman)))
  
}


options(digits = 2)
get_LRs(input=dat.450_540.genotype, container = "CBS Straw Green")
get_LRs(input=dat.450_540.genotype, container = "CBS Straw White")

get_LRs(input=dat.450_570.genotype, container = "CBS Straw Green")
get_LRs(input=dat.450_570.genotype, container = "CBS Straw White")


# please run linear regression between Y=concentration of cleaved Par4 and X=log(0.5*(30K concentration)) 
# as well as pearson and spearman correlation coefficients 
# for all samples, genotype=0, genotype=1, genotype=1 or 2 samples respectively after removing that one sample.

get_LRs <- function(input, container) {
  #----------------------------------------#
  # Y=concentration of cleaved Par4 and X=log(0.5*(30K concentration)) 
  coef.pearson <- c(); coef.spearman <- c()
  ## all samples
  y <- input$cleaved_par4[input$Container == container]; x <- log(0.5*input$X30K[input$Container == container])
  LR.sum <- summary(lm(y ~ x))
  LR.out.all <- LR.sum$coefficients
  coef.pearson.all <- cor(x, y, method = c("pearson"))
  coef.spearman.all <- cor(x, y, method = c("spearman"))
  
  ## genotype=0
  y <- input$cleaved_par4[input$rs6494889==0 & input$Container == container]; 
  x <- log(0.5*input$X30K)[input$rs6494889==0 & input$Container == container] 
  LR.sum <- summary(lm(y ~ x))
  LR.out.0 <- LR.sum$coefficients
  coef.pearson.0 <- cor(x, y, method = c("pearson"))
  coef.spearman.0 <- cor(x, y, method = c("spearman"))
  
  ## genotype=1
  y <- input$cleaved_par4[input$rs6494889==1 & input$Container == container]; 
  x <- log(0.5*input$X30K)[input$rs6494889==1 & input$Container == container] 
  LR.sum <- summary(lm(y ~ x))
  LR.out.1 <- LR.sum$coefficients
  coef.pearson.1 <- cor(x, y, method = c("pearson"))
  coef.spearman.1 <- cor(x, y, method = c("spearman"))
  
  ## genotype=1+2
  y <- input$cleaved_par4[input$rs6494889!=0 & input$Container == container]; 
  x <- log(0.5*input$X30K)[input$rs6494889!=0 & input$Container == container] 
  LR.sum <- summary(lm(y ~ x))
  LR.out.12 <- LR.sum$coefficients
  coef.pearson.12 <- cor(x, y, method = c("pearson"))
  coef.spearman.12 <- cor(x, y, method = c("spearman"))
  
  # output
  LR <- rbind(LR.out.all, LR.out.0, LR.out.1, LR.out.12)
  coef.pearson <- c(coef.pearson.all, coef.pearson.0, coef.pearson.1, coef.pearson.12)
  coef.spearman <- c(coef.spearman.all,  coef.spearman.0, coef.spearman.1, coef.spearman.12)
  
  return(list("LR"=LR, 
              "coef.pearson"=data.frame("name"=c("all", "0", "1", "1+2"), "coef"=coef.pearson),
              "coef.spearman"=data.frame("name"=c("all", "0", "1", "1+2"), "coef"=coef.spearman)))
  
}

options(digits = 2)
get_LRs(input=dat.450_540.genotype, container = "CBS Straw Green")
get_LRs(input=dat.450_540.genotype, container = "CBS Straw White")

get_LRs(input=dat.450_570.genotype, container = "CBS Straw Green")
get_LRs(input=dat.450_570.genotype, container = "CBS Straw White")



# please also run regression for Y=X+container*genotype, where Y=concentration of cleaved Par4 and X=0.5*(30K concentration).
get_LRs <- function(input) {
  input$y <- input$cleaved_par4
  input$x <- 0.5*input$X30K
  input$rs6494889 <- as.factor(input$rs6494889)
  input = within(input, Container <- relevel(as.factor(Container), ref = "CBS Straw White"))
  
  LR.sum <- summary(lm(y ~ x + Container*rs6494889, data=input))
  return(LR.sum)
}
get_LRs(input=dat.450_540.genotype)
get_LRs(input=dat.450_570.genotype)

# when you run  regression for Y=X+container*genotype, where Y=concentration of cleaved Par4 and X=0.5*(30K concentration), 
# please treat genotype as a numeric variable
get_LRs <- function(input) {
  input$y <- input$cleaved_par4
  input$x <- 0.5*input$X30K
  input$rs6494889 <- as.numeric(input$rs6494889)
  input = within(input, Container <- relevel(as.factor(Container), ref = "CBS Straw White"))
  
  LR.sum <- summary(lm(y ~ x + Container*rs6494889, data=input))
  return(LR.sum)
}
get_LRs(input=dat.450_540.genotype)
get_LRs(input=dat.450_570.genotype)

#----------------------------------------------------------------------------------------------

input=dat.450_540.genotype
# y <- input$cleaved_par4[input$rs6494889!=0]; x <- 0.5*input$X30K[input$rs6494889!=0] 
# group <- input$Container[input$rs6494889!=0]
y <- input$cleaved_par4[input$rs6494889==0]; x <- 0.5*input$X30K[input$rs6494889==0] 
group <- input$Container[input$rs6494889==0]
df <- data.frame("x"=x, "y"=y, "group"=group)
setwd("/user/xinyugao/ELISA_project/")
pdf(paste0("cleavedPar4_vs_half30K_450_540nm.genotype_1or2.pdf"))
plot(x, y,
     pch = 19,
     col = factor(group), 
     ylab="cleaved_par4 concentration", xlab = "0.5 * 30K concentration")
LR.sum <- summary(lm(y ~ x))
abline(a = LR.sum$coefficients[1,1], b = LR.sum$coefficients[2,1], col="blue")
LR.green <- summary(lm(y ~ x, data = df, subset = group=="CBS Straw Green"))
abline(a = LR.green$coefficients[1,1], b = LR.green$coefficients[2,1], col="black", lty=2)
LR.white <- summary(lm(y ~ x, data = df, subset = group=="CBS Straw White"))
abline(a = LR.white$coefficients[1,1], b = LR.white$coefficients[2,1], col="red", lty=2)
# Legend
legend("topright",
       legend = levels(factor(group)),
       pch = 19,
       col = factor(levels(factor(group))))
dev.off()


input=dat.450_540.genotype
y <- input$original[input$rs6494889!=0]; x <- 0.5*input$X30K[input$rs6494889!=0]
group <- input$Container[input$rs6494889!=0]
# tmp <- input %>% select(original, X30K, Container) %>% filter(input$rs6494889!=0)
# y <- input$original[input$rs6494889==0]; x <- 0.5*input$X30K[input$rs6494889==0]
# group <- input$Container[input$rs6494889==0]
df <- data.frame("x"=x, "y"=y, "group"=group)
setwd("/user/xinyugao/ELISA_project/")
pdf(paste0("original_vs_half30K_450_540nm.genotype_1or2.pdf"))
# y <- input$original; x <- 0.5*input$X30K
plot(x, y,
     pch = 19,
     col = factor(group),
     ylab="original concentration", xlab = "0.5 * 30K concentration")
LR.sum <- summary(lm(y ~ x))
abline(a = LR.sum$coefficients[1,1], b = LR.sum$coefficients[2,1], col="blue")
LR.green <- summary(lm(y ~ x, data = df, subset = group=="CBS Straw Green"))
abline(a = LR.green$coefficients[1,1], b = LR.green$coefficients[2,1], col="black", lty=2)
LR.white <- summary(lm(y ~ x, data = df, subset = group=="CBS Straw White"))
abline(a = LR.white$coefficients[1,1], b = LR.white$coefficients[2,1], col="red", lty=2)
# Legend
legend("topright",
       legend = levels(factor(group)),
       pch = 19,
       col = factor(levels(factor(group))))
dev.off()



# please run linear regression between Y=concentration of cleaved Par4 and X=4^(0.5*(30K concentration)) 
# as well as pearson and spearman correlation coefficients for 
# all samples, genotype=0, genotype=1, genotype=1 or 2 samples respectively after removing that one sample.


get_LRs <- function(input, coef.val=1.064) {
  #----------------------------------------#
  # Y=concentration of cleaved Par4 and X=4^(0.5*(30K concentration))  y = coef* x^4
  coef.pearson <- c(); coef.spearman <- c()
  ## all samples
  # y <- input$cleaved_par4; x <- 4^(0.5*input$X30K)
  # y <- input$cleaved_par4; x <- (0.5*input$X30K)^4
  y <- input$cleaved_par4; x <-  coef.val^(0.5*input$X30K)

  LR.sum <- summary(lm(y ~ x))
  LR.out.all <- LR.sum$coefficients
  coef.pearson.all <- cor(x, y, method = c("pearson"))
  coef.spearman.all <- cor(x, y, method = c("spearman"))

  ## genotype=0
  # y <- input$cleaved_par4[input$rs6494889==0]; x <- 4^(0.5*input$X30K)[input$rs6494889==0]
  y <- input$cleaved_par4[input$rs6494889==0]; x <- (coef.val^(0.5*input$X30K))[input$rs6494889==0]
  
  # y ~ coef * b_0 ^(0.5*30K) b_0 using LS

  LR.sum <- summary(lm(y ~ x))
  LR.out.0 <- LR.sum$coefficients
  coef.pearson.0 <- cor(x, y, method = c("pearson"))
  coef.spearman.0 <- cor(x, y, method = c("spearman"))

  ## genotype=1
  # y <- input$cleaved_par4[input$rs6494889==1]; x <- 4^(0.5*input$X30K)[input$rs6494889==1]
  y <- input$cleaved_par4[input$rs6494889==1]; x <- (coef.val^(0.5*input$X30K))[input$rs6494889==1]
  LR.sum <- summary(lm(y ~ x))
  LR.out.1 <- LR.sum$coefficients
  coef.pearson.1 <- cor(x, y, method = c("pearson"))
  coef.spearman.1 <- cor(x, y, method = c("spearman"))

  ## genotype=1+2
  # y <- input$cleaved_par4[input$rs6494889!=0]; x <- 4^(0.5*input$X30K)[input$rs6494889!=0]
  y <- input$cleaved_par4[input$rs6494889!=0]; x <- (coef.val^(0.5*input$X30K))[input$rs6494889!=0]
  
  LR.sum <- summary(lm(y ~ x))
  LR.out.12 <- LR.sum$coefficients
  coef.pearson.12 <- cor(x, y, method = c("pearson"))
  coef.spearman.12 <- cor(x, y, method = c("spearman"))

  # output
  LR <- rbind(LR.out.all, LR.out.0, LR.out.1, LR.out.12)
  coef.pearson <- c(coef.pearson.all, coef.pearson.0, coef.pearson.1, coef.pearson.12)
  coef.spearman <- c(coef.spearman.all,  coef.spearman.0, coef.spearman.1, coef.spearman.12)

  return(list("LR"=LR,
              "coef.pearson"=data.frame("name"=c("all", "0", "1", "1+2"), "coef"=coef.pearson),
              "coef.spearman"=data.frame("name"=c("all", "0", "1", "1+2"), "coef"=coef.spearman)))

}

options(digits = 2)
get_LRs(input=dat.450_540.genotype, coef.val = 1.064)
get_LRs(input=dat.450_570.genotype)

input=dat.450_570.genotype
y <- input$cleaved_par4; x <- (0.5*input$X30K)
# y <- input$cleaved_par4[input$rs6494889!=0]; x <- (0.5*input$X30K)[input$rs6494889!=0]
my.formula <- formula(y~b^x)
model.fit <- nls(my.formula, start = list(b=1))
sum.fit <- summary(model.fit)
options(digits = 5)
sum.fit$coefficients
y.pred <- predict(model.fit)
plot(x, y)
points(x, y.pred,col="red",lwd=1, pch=3)




# (original, 30K, cleaved par4)
# Y=Par-4 level ~ DxAge + BMI2 + as.factor(ER) + as.factor(PR) + as.factor(HER2) + as.factor(Hormonal.Therapy) +
#   as.factor(Radiation.Therapy) + as.factor(Grade_Description) + as.factor(StageNEW) + as.factor(Surgery.Type) + rs6494889 + (COLLECTIONDATE-DateDx)

head(PT.genotype.rm)
colSums(is.na(PT.genotype.rm))
pheno <- PT.genotype.rm

pheno = within(pheno, ER <- relevel(as.factor(ER), ref = "2"))
pheno = within(pheno, PR <- relevel(as.factor(PR), ref = "2"))
pheno = within(pheno, HER2 <- relevel(as.factor(HER2), ref = "2"))
pheno$DxAge = as.numeric(as.character(pheno$DxAge))
#pheno = within(pheno, Surgery.Type.clean <- relevel(as.factor(Surgery.Type.clean), ref = "None"))
pheno = within(pheno, Surgery.Type <- relevel(as.factor(Surgery.Type), ref = "Lumpectomy"))
pheno$rs6494889 = as.numeric(as.character(pheno$rs6494889))
pheno$datediff <- difftime(as.Date(pheno$COLLECTIONDATE, "%Y-%m-%d"), as.Date(pheno$DateDx, format = "%m/%d/%Y"), units = "days")


dat.450_540.merged <- dat.450_540.genotype %>% select(-rs6494889) 
dat.450_570.merged <- dat.450_570.genotype %>% select(-rs6494889)

pheno.450_540 <- inner_join(pheno, dat.450_540.merged, by=c("SubjectID"="Subjectid"))
pheno.450_570 <- inner_join(pheno, dat.450_570.merged, by=c("SubjectID"="Subjectid"))
pheno.450_540$log_original <- log(pheno.450_540$original)
pheno.450_540$log_30K <- log(pheno.450_540$X30K)
pheno.450_570$log_original <- log(pheno.450_570$original)
pheno.450_570$log_30K <- log(pheno.450_570$X30K)


res = lm(original ~ DxAge + BMI2 + as.factor(ER) + as.factor(PR) + as.factor(HER2) + as.factor(Hormonal.Therapy) +
           as.factor(Radiation.Therapy) + as.factor(Grade_Description) + as.factor(StageNEW) + 
           as.factor(Surgery.Type) + rs6494889 + datediff, data=pheno.450_570) 
summary(res)$coefficients


#####################################################################################################################
# please generate another variable called HR. HR is 1 if either ER or PR is 1; HR is missing if either ER or PR is missing; otherwise HR is 0.

head(PT.genotype.rm)
colSums(is.na(PT.genotype.rm))
pheno <- PT.genotype.rm

pheno = within(pheno, ER <- relevel(as.factor(ER), ref = "2"))
pheno = within(pheno, PR <- relevel(as.factor(PR), ref = "2"))
pheno = within(pheno, HER2 <- relevel(as.factor(HER2), ref = "2"))
pheno$DxAge = as.numeric(as.character(pheno$DxAge))
#pheno = within(pheno, Surgery.Type.clean <- relevel(as.factor(Surgery.Type.clean), ref = "None"))
pheno = within(pheno, Surgery.Type <- relevel(as.factor(Surgery.Type), ref = "Lumpectomy"))
pheno$rs6494889 = as.numeric(as.character(pheno$rs6494889))
pheno$datediff <- difftime(as.Date(pheno$COLLECTIONDATE, "%Y-%m-%d"), as.Date(pheno$DateDx, format = "%m/%d/%Y"), units = "days")


pheno$HR <- NA
for (i in 1:dim(pheno)[1]) {
  if (is.na(pheno$ER[i]) | is.na(pheno$PR[i])) {
    pheno$HR <- NA
  } else if (pheno$ER[i] == 1 | pheno$PR[i] == 1) {
    pheno$HR[i] <- 1
  } else {
    pheno$HR[i] <- 0
  }
}

# please generate 2x2 table for HR and Hormonal.Therapy and run fisher exact test
table(pheno$HR)
# 0  1 
# 14 24 
table(pheno$Hormonal.Therapy)
# 0  1 
# 13 25 
table(pheno$HR, pheno$Hormonal.Therapy)
#    0  1
# 0 12  2
# 1  1 23
df <- data.frame("Hormonal.Therapy=0" = c(12, 1), "Hormonal.Therapy=1" = c(2, 23), row.names = c("HR=0", "HR=1"))
df
mosaicplot(df, color = TRUE, "mosaicplot")  
library(stats)
fisher.test(df)
# Fisher's Exact Test for Count Data
# 
# data:  df
# p-value = 4.059e-07
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#     9.157196 5802.724304
# sample estimates:
# odds ratio 
#   102.6059 


dat.450_540.merged <- dat.450_540.genotype %>% select(-rs6494889) 
dat.450_570.merged <- dat.450_570.genotype %>% select(-rs6494889)

pheno.450_540 <- inner_join(pheno, dat.450_540.merged, by=c("SubjectID"="Subjectid"))
pheno.450_570 <- inner_join(pheno, dat.450_570.merged, by=c("SubjectID"="Subjectid"))
pheno.450_540$log_original <- log(pheno.450_540$original)
pheno.450_540$log_30K <- log(pheno.450_540$X30K)
pheno.450_570$log_original <- log(pheno.450_570$original)
pheno.450_570$log_30K <- log(pheno.450_570$X30K)

options(digits = 3)
res = lm(original ~ DxAge + BMI2 + as.factor(HR) + as.factor(Hormonal.Therapy) + as.factor(HER2) + 
           as.factor(Radiation.Therapy) + as.factor(Grade_Description) + as.factor(StageNEW) +  
           as.factor(Surgery.Type) + rs6494889 + datediff, data=pheno.450_570) 
summary(res)$coefficients


# alculate COLLECTIONDATE-Hormonal.Date and also include this variable in the regression model
# add to this regression model Y=Par-4 level ~ DxAge + BMI2 + as.factor(ER) + as.factor(PR) + as.factor(HER2) + as.factor(Hormonal.Therapy) +
#   as.factor(Radiation.Therapy) + as.factor(Grade_Description) + as.factor(StageNEW) + as.factor(Surgery.Type) + rs6494889 + (COLLECTIONDATE-DateDx)
head(PT.genotype.rm)
colSums(is.na(PT.genotype.rm))
pheno <- 

pheno = within(pheno, ER <- relevel(as.factor(ER), ref = "2"))
pheno = within(pheno, PR <- relevel(as.factor(PR), ref = "2"))
pheno = within(pheno, HER2 <- relevel(as.factor(HER2), ref = "2"))
pheno$DxAge = as.numeric(as.character(pheno$DxAge))
#pheno = within(pheno, Surgery.Type.clean <- relevel(as.factor(Surgery.Type.clean), ref = "None"))
pheno = within(pheno, Surgery.Type <- relevel(as.factor(Surgery.Type), ref = "Lumpectomy"))
pheno$rs6494889 = as.numeric(as.character(pheno$rs6494889))
pheno$datediff <- as.numeric(difftime(as.Date(pheno$COLLECTIONDATE, "%Y-%m-%d"), as.Date(pheno$DateDx, format = "%m/%d/%Y"), units = "days"))
pheno$datediff2 <- as.numeric(difftime(as.Date(pheno$Hormonal.Date, format = "%m/%d/%Y"), as.Date(pheno$COLLECTIONDATE, "%Y-%m-%d"), units = "days"))

dat.450_540.merged <- dat.450_540.genotype %>% select(-rs6494889) 
dat.450_570.merged <- dat.450_570.genotype %>% select(-rs6494889)

pheno.450_540 <- inner_join(pheno, dat.450_540.merged, by=c("SubjectID"="Subjectid"))
pheno.450_570 <- inner_join(pheno, dat.450_570.merged, by=c("SubjectID"="Subjectid"))
pheno.450_540$log_original <- log(pheno.450_540$original)
pheno.450_540$log_30K <- log(pheno.450_540$X30K)
pheno.450_570$log_original <- log(pheno.450_570$original)
pheno.450_570$log_30K <- log(pheno.450_570$X30K)

pheno[pheno$HR == 0 & pheno$Hormonal.Therapy == 1, ] %>% select(c("SubjectID", "COLLECTIONBARCODE", "CollectionID", "HR", "Hormonal.Therapy","HER2"))
# SubjectID     COLLECTIONBARCODE CollectionID HR   Hormonal.Therapy  HER2
# PT-00286063         132131452   Cl-00224002  0                1      2
# PT-00192391          72281664   Cl-00039089  0                1      1
pheno[pheno$HR == 1 & pheno$Hormonal.Therapy == 0, ] %>% select(c("SubjectID", "COLLECTIONBARCODE", "CollectionID", "HR", "Hormonal.Therapy","HER2"))
# SubjectID    COLLECTIONBARCODE CollectionID   HR   Hormonal.Therapy  HER2
# PT-00216317          93381175  Cl-00046326     1           0           1


# tmp <- pheno %>% select(DxAge, BMI2 ,ER, PR, HER2, Hormonal.Therapy,
#                         Radiation.Therapy, Grade_Description, StageNEW, Surgery.Type, rs6494889 , datediff , datediff2)
# tmp2 <- na.omit(tmp) 
# dim(tmp2)

options(digits = 3)
res = lm(log_30K~ DxAge + BMI2 + as.factor(ER) + as.factor(PR) + as.factor(HER2) + 
           as.factor(Radiation.Therapy) + as.factor(Grade_Description) + as.factor(StageNEW) + 
           as.factor(Surgery.Type) + rs6494889 + datediff + datediff2, data=pheno.450_570) 

summary(res)$coefficients



# please add two additional covariates ("Container" from DistSamplesList_nowaka_Sep-16-2022-11-25-05.xls and today-COLLECTIONDATE) 
# to this regression model Y=Par-4 level ~ DxAge + BMI2 + as.factor(ER) + as.factor(PR) + as.factor(HER2) + as.factor(Hormonal.Therapy) +
#   as.factor(Radiation.Therapy) + as.factor(Grade_Description) + as.factor(StageNEW) + as.factor(Surgery.Type) + rs6494889 + (COLLECTIONDATE-DateDx)


head(PT.genotype.rm)
colSums(is.na(PT.genotype.rm))
pheno <- PT.genotype.rm

pheno = within(pheno, ER <- relevel(as.factor(ER), ref = "2"))
pheno = within(pheno, PR <- relevel(as.factor(PR), ref = "2"))
pheno = within(pheno, HER2 <- relevel(as.factor(HER2), ref = "2"))

pheno$DxAge = as.numeric(as.character(pheno$DxAge))
#pheno = within(pheno, Surgery.Type.clean <- relevel(as.factor(Surgery.Type.clean), ref = "None"))
pheno = within(pheno, Surgery.Type <- relevel(as.factor(Surgery.Type), ref = "Lumpectomy"))
pheno$rs6494889 = as.numeric(as.character(pheno$rs6494889))
pheno$datediff <- difftime(as.Date(pheno$COLLECTIONDATE, "%Y-%m-%d"), as.Date(pheno$DateDx, format = "%m/%d/%Y"), units = "days")
pheno$datediff3 <- difftime("2022-09-30", as.Date(pheno$COLLECTIONDATE, "%Y-%m-%d"), units = "days")

pheno$HR <- NA
for (i in 1:dim(pheno)[1]) {
  if (is.na(pheno$ER[i]) | is.na(pheno$PR[i])) {
    pheno$HR <- NA
  } else if (pheno$ER[i] == 1 | pheno$PR[i] == 1) {
    pheno$HR[i] <- 1
  } else {
    pheno$HR[i] <- 0
  }
}

dat.450_540.merged <- dat.450_540.genotype %>% select(-rs6494889) 
dat.450_570.merged <- dat.450_570.genotype %>% select(-rs6494889)

pheno.450_540 <- inner_join(pheno, dat.450_540.merged, by=c("SubjectID"="Subjectid"))
pheno.450_570 <- inner_join(pheno, dat.450_570.merged, by=c("SubjectID"="Subjectid"))
pheno.450_540$log_original <- log(pheno.450_540$original)
pheno.450_540$log_30K <- log(pheno.450_540$X30K)
pheno.450_570$log_original <- log(pheno.450_570$original)
pheno.450_570$log_30K <- log(pheno.450_570$X30K)


# pheno.450_540 = within(pheno.450_540, Container <- relevel(as.factor(Container), ref = "CBS Straw White"))
# pheno.450_570 = within(pheno.450_570, Container <- relevel(as.factor(Container), ref = "CBS Straw White"))
# options(digits = 3)
# res = lm(original ~ DxAge + BMI2 + as.factor(ER) + as.factor(PR) + as.factor(HER2) + as.factor(Hormonal.Therapy) +
#            as.factor(Radiation.Therapy) + as.factor(Grade_Description) + as.factor(StageNEW) + 
#            as.factor(Surgery.Type) + rs6494889 + datediff + datediff3 + as.factor(Container), 
#          data=pheno.450_570) 
# summary(res)$coefficients


options(digits = 3)
res = lm(original ~ DxAge + BMI2 + as.factor(ER) + as.factor(PR) + as.factor(HER2) + as.factor(Hormonal.Therapy) +
           as.factor(Radiation.Therapy) + as.factor(Grade_Description) + as.factor(StageNEW) + 
           as.factor(Surgery.Type) + rs6494889 + datediff + datediff3, 
          data=pheno.450_540, subset = Container == "CBS Straw Green") 
summary(res)$coefficients


# tt <- pheno.450_540[pheno.450_540$Container == "CBS Straw Green", c("DxAge", "BMI2", "ER", "PR", "HER2",  "Hormonal.Therapy","Radiation.Therapy","Grade_Description",
#   "StageNEW", "Surgery.Type", "rs6494889", "datediff", "datediff3")]
# tt2 <- pheno.450_540[pheno.450_540$Container == "CBS Straw White", c("DxAge", "BMI2", "ER", "PR", "HER2",  "Hormonal.Therapy","Radiation.Therapy","Grade_Description",
#                                                                      "StageNEW", "Surgery.Type", "rs6494889", "datediff", "datediff3")]

# options(digits = 3)
# res = lm(log_30K ~ DxAge + BMI2 + as.factor(ER) + as.factor(PR) + as.factor(HER2) + as.factor(Hormonal.Therapy) +
#            as.factor(Radiation.Therapy) + as.factor(Grade_Description) + as.factor(StageNEW) + 
#            as.factor(Surgery.Type) + rs6494889 + datediff + datediff3, 
#          data=pheno.450_570, subset = Container == "CBS Straw White") 
# summary(res)$coefficients


options(digits = 3)
res = lm(log_30K ~ DxAge + BMI2 + as.factor(HR) + as.factor(Hormonal.Therapy) + as.factor(HER2) + 
           as.factor(Radiation.Therapy) + as.factor(Grade_Description) + as.factor(StageNEW) + 
           as.factor(Surgery.Type) + rs6494889 + datediff + datediff3, 
         data=pheno.450_570, subset = Container == "CBS Straw White") 
summary(res)$coefficients


##########################################################################################################################
# please run model 5.2:  Par-4 level ~DxAge + BMI2 + as.factor(ER) + as.factor(PR) + as.factor(HER2) + as.factor(Hormonal.Therapy)
# +as.factor(Radiation.Therapy) + as.factor(Grade_Description) + as.factor(StageNEW) + as.factor(Surgery.Type) + datediff + datediff3 + Container*rs6494889
# treat rs6494889 as a numeric variable


head(PT.genotype.rm)
colSums(is.na(PT.genotype.rm))
pheno <- PT.genotype.rm

pheno = within(pheno, ER <- relevel(as.factor(ER), ref = "2"))
pheno = within(pheno, PR <- relevel(as.factor(PR), ref = "2"))
pheno = within(pheno, HER2 <- relevel(as.factor(HER2), ref = "2"))

pheno$DxAge = as.numeric(as.character(pheno$DxAge))
#pheno = within(pheno, Surgery.Type.clean <- relevel(as.factor(Surgery.Type.clean), ref = "None"))
pheno = within(pheno, Surgery.Type <- relevel(as.factor(Surgery.Type), ref = "Lumpectomy"))
pheno$rs6494889 = as.numeric(as.character(pheno$rs6494889))
pheno$datediff <- difftime(as.Date(pheno$COLLECTIONDATE, "%Y-%m-%d"), as.Date(pheno$DateDx, format = "%m/%d/%Y"), units = "days")
pheno$datediff3 <- difftime("2022-09-30", as.Date(pheno$COLLECTIONDATE, "%Y-%m-%d"), units = "days")

dat.450_540.merged <- dat.450_540.genotype %>% select(-rs6494889) 
dat.450_570.merged <- dat.450_570.genotype %>% select(-rs6494889)

pheno.450_540 <- inner_join(pheno, dat.450_540.merged, by=c("SubjectID"="Subjectid"))
pheno.450_570 <- inner_join(pheno, dat.450_570.merged, by=c("SubjectID"="Subjectid"))
pheno.450_540$log_original <- log(pheno.450_540$original)
pheno.450_540$log_30K <- log(pheno.450_540$X30K)
pheno.450_570$log_original <- log(pheno.450_570$original)
pheno.450_570$log_30K <- log(pheno.450_570$X30K)

pheno.450_540 = within(pheno.450_540, Container <- relevel(as.factor(Container), ref = "CBS Straw White"))
pheno.450_570 = within(pheno.450_570, Container <- relevel(as.factor(Container), ref = "CBS Straw White"))

options(digits = 3)
res = lm(log_30K ~ DxAge + BMI2 + as.factor(ER) + as.factor(PR) + as.factor(HER2) + as.factor(Hormonal.Therapy) + 
           as.factor(Radiation.Therapy) + as.factor(Grade_Description) + as.factor(StageNEW) + 
           as.factor(Surgery.Type) + Container*rs6494889 + datediff + datediff3, data=pheno.450_570) 
summary(res)$coefficients



###########################################################################################################################
# 10/19 TASK: model 7: Par-4 level ~DxAge + BMI2 + as.factor(ER) + as.factor(PR) + as.factor(HER2) + as.factor(Grade_Description) + 
# as.factor(StageNEW) + as.factor(Surgery.Type) + datediff + datediff3 + Container+rs6494889?    treat rs6494889 as a numeric variable.
# please run model 7.1: Par-4 level ~DxAge + BMI2 + as.factor(ER) + as.factor(PR) + as.factor(HER2) + as.factor(Grade_Description) + 
# as.factor(StageNEW) + as.factor(Surgery.Type) + datediff + datediff3 + Container*rs6494889? treat rs6494889 as a numeric variable.

head(PT.genotype.rm)
colSums(is.na(PT.genotype.rm))
pheno <- PT.genotype.rm

PT.genotype.rm$COLLECTIONDATE
PT.genotype.rm$Surgery.Date

# tt <- PT.genotype.rm %>% select(SubjectID, COLLECTIONDATE, Surgery.Date)
# tt$COLLECTIONDATE <- as.Date(tt$COLLECTIONDATE, "%Y-%m-%d")
# tt$Surgery.Date <- as.Date(tt$Surgery.Date, "%m/%d/%Y")
# tt$collectionBeforeSurgery <- tt$COLLECTIONDATE < tt$Surgery.Date
# table(tt$collectionBeforeSurgery)

pheno = within(pheno, ER <- relevel(as.factor(ER), ref = "2"))
pheno = within(pheno, PR <- relevel(as.factor(PR), ref = "2"))
pheno = within(pheno, HER2 <- relevel(as.factor(HER2), ref = "2"))

pheno$DxAge = as.numeric(as.character(pheno$DxAge))
#pheno = within(pheno, Surgery.Type.clean <- relevel(as.factor(Surgery.Type.clean), ref = "None"))
pheno = within(pheno, Surgery.Type <- relevel(as.factor(Surgery.Type), ref = "Lumpectomy"))
pheno$rs6494889 = as.numeric(as.character(pheno$rs6494889))
pheno$datediff <- difftime(as.Date(pheno$COLLECTIONDATE, "%Y-%m-%d"), as.Date(pheno$DateDx, format = "%m/%d/%Y"), units = "days")
pheno$datediff3 <- difftime("2022-09-30", as.Date(pheno$COLLECTIONDATE, "%Y-%m-%d"), units = "days")

dat.450_540.merged <- dat.450_540.genotype %>% select(-rs6494889) 
dat.450_570.merged <- dat.450_570.genotype %>% select(-rs6494889)

pheno.450_540 <- inner_join(pheno, dat.450_540.merged, by=c("SubjectID"="Subjectid"))
pheno.450_570 <- inner_join(pheno, dat.450_570.merged, by=c("SubjectID"="Subjectid"))
pheno.450_540$log_original <- log(pheno.450_540$original)
pheno.450_540$log_30K <- log(pheno.450_540$X30K)
pheno.450_570$log_original <- log(pheno.450_570$original)
pheno.450_570$log_30K <- log(pheno.450_570$X30K)

pheno.450_540 = within(pheno.450_540, Container <- relevel(as.factor(Container), ref = "CBS Straw White"))
pheno.450_570 = within(pheno.450_570, Container <- relevel(as.factor(Container), ref = "CBS Straw White"))

# options(digits = 3)
# this_data <- pheno.450_570
# # M7: no interaction
# res = lm(log_30K ~ DxAge + BMI2 + as.factor(ER) + as.factor(PR) + as.factor(HER2) +  as.factor(Grade_Description) + as.factor(StageNEW) +
#             Container+rs6494889 + datediff , data=this_data)
# summary(res)$coefficients
# 
# # M7.1: with interaction
# res = lm(log_30K ~ DxAge + BMI2 + as.factor(ER) + as.factor(PR) + as.factor(HER2) +  as.factor(Grade_Description) + as.factor(StageNEW) +
#            Container*rs6494889 + datediff , data=this_data)
# summary(res)$coefficients

# M7.2
options(digits = 3)
this_data <- pheno.450_570
res = lm(cleaved_par4 ~ DxAge + BMI2 + as.factor(ER) + as.factor(PR) + as.factor(HER2) +  as.factor(Grade_Description) + as.factor(StageNEW) +
           rs6494889 + datediff , data=this_data, subset = Container == "CBS Straw White")
summary(res)$coefficients

# M7.3
res = lm(original ~ DxAge + BMI2 + as.factor(ER) + as.factor(PR) + as.factor(HER2) +  as.factor(Grade_Description) + as.factor(StageNEW) +
           rs6494889 + datediff , data=this_data, subset = Container == "CBS Straw Green")
summary(res)$coefficients

# tt <- pheno.450_540[pheno.450_540$Container == "CBS Straw Green", c("DxAge", "BMI2", "ER", "PR", "HER2","Grade_Description",
#   "StageNEW",  "rs6494889", "datediff", "datediff3")]
# tt2 <- pheno.450_540[pheno.450_540$Container == "CBS Straw White", c("DxAge", "BMI2", "ER", "PR", "HER2","Grade_Description",
#                                                                      "StageNEW",  "rs6494889", "datediff", "datediff3")]








###########################################################################################################################
library(readxl)

filepath.4th.concentration.540nm <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/inputs/4th/4thELISA_450nm-540nm"
filepath.4th.concentration.570nm <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/inputs/4th/4thELISA_450nm-570nm"
filepath.4th.concentration.540nm.30K <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/inputs/4th/4thELISA_450nm-540nm.30K"
filepath.4th.concentration.570nm.30K <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/inputs/4th/4thELISA_450nm-570nm.30K"

filepath.4th.sample.info <- "/projects/rpci/qzhu/xinyu/mywork/ELISA_par4/inputs/4th/4th_ELISA_samples_results.xlsx"# sample info

sample.4th.info <- read_excel(path = filepath.4th.sample.info, sheet = "Sheet1")

# The "Plasma Type" column in /projects/rpci/qzhu/users-qzhu/grant/UACA_R01/ELISA/Yao Available Blood 11.1.21_wTube.xlsx corresponds to container information.
# White straw from edta vacutainer
# Green straw from heparin vacutainer


container.5th.info <- readxl::read_excel("/projects/rpci/qzhu/users-qzhu/grant/UACA_R01/ELISA/Yao Available Blood 11.1.21_wTube.xlsx") # collection_alias



concen.4th.540.original <- read.table(file = filepath.4th.concentration.540nm, header = TRUE, colClasses = c("numeric", "numeric"))
concen.4th.540.30K <- read.table(file = filepath.4th.concentration.540nm.30K, header = TRUE, colClasses = c("numeric", "numeric"))
concen.4th.570.original <- read.table(file = filepath.4th.concentration.570nm, header = TRUE, colClasses = c("numeric", "numeric"))
concen.4th.570.30K <- read.table(file = filepath.4th.concentration.570nm.30K, header = TRUE, colClasses = c("numeric", "numeric"))

concen.4th.540.all <- inner_join(concen.4th.540.original, concen.4th.540.30K, by="ID")
colnames(concen.4th.540.all) <- c("ID" , "original", "X30K")
concen.4th.570.all <- inner_join(concen.4th.570.original, concen.4th.570.30K, by="ID")
colnames(concen.4th.570.all) <- c("ID" , "original", "X30K")

tmp.540 <- inner_join(sample.4th.info, concen.4th.540.all, by=c("BARCODE"="ID"))
tmp.570 <- inner_join(sample.4th.info, concen.4th.570.all, by=c("BARCODE"="ID"))

dat_reg.540 <- inner_join(tmp.540, container.5th.info, by=c("SubjectID", "COLLECTIONID" ))
dat_reg.570 <- inner_join(tmp.570, container.5th.info, by=c("SubjectID", "COLLECTIONID" ))
dat_reg.540$cleaved_par4 <- dat_reg.540$original - 0.5*dat_reg.540$X30K
dat_reg.570$cleaved_par4 <- dat_reg.570$original - 0.5*dat_reg.570$X30K

xy <- dat_reg.540 %>% select(cleaved_par4, `Plasma Type`)


get_LR <- function(input, type="Heparin") {
  coef.pearson <- c(); coef.spearman <- c()
  ## all samples 
  # original vs 0.5*30K
  y <- input$original[input$`Plasma Type` == type]; x <- 0.5 * input$X30K[input$`Plasma Type` == type]
  LR.sum <- summary(lm(y ~ x))
  LR.out.all.mod1 <- LR.sum$coefficients
  coef.pearson.all.mod1 <- cor(x, y, method = c("pearson"), use="complete.obs")
  coef.spearman.all.mod1 <- cor(x, y, method = c("spearman"), use="complete.obs")
  
  ## all samples 
  # cleaved_par4 vs 0.5*30K
  y <- input$cleaved_par4[input$`Plasma Type` == type]; x <- 0.5 * input$X30K[input$`Plasma Type` == type]
  LR.sum <- summary(lm(y ~ x))
  LR.out.all.mod2 <- LR.sum$coefficients
  coef.pearson.all.mod2 <- cor(x, y, method = c("pearson"), use="complete.obs")
  coef.spearman.all.mod2 <- cor(x, y, method = c("spearman"), use="complete.obs")
  
  
  # output
  LR <- rbind(LR.out.all.mod1, LR.out.all.mod2)
  coef.pearson <- c(coef.pearson.all.mod1, coef.pearson.all.mod2)
  coef.spearman <- c(coef.spearman.all.mod1, coef.spearman.all.mod2)
  
  return(list("LR"=LR,
              "coef.pearson"=data.frame("name"=c("mod1", "mod2"), "coef"=coef.pearson),
              "coef.spearman"=data.frame("name"=c("mod1", "mod2"), "coef"=coef.spearman)))
}

options(digits = 3)
get_LR(input = dat_reg.540)
get_LR(input = dat_reg.570)

input = dat_reg.540
# y <- input$cleaved_par4[input$`Plasma Type` == "Heparin"]; x <- 0.5 * input$X30K[input$`Plasma Type` == "Heparin"]
# plot(x, y, col="black", pch=19, xlab = "0.5*30K", ylab="cleaved_par4")
y <- input$original[input$`Plasma Type` == "Heparin"]; x <- 0.5 * input$X30K[input$`Plasma Type` == "Heparin"]
plot(x, y, col="black", pch=19, xlab = "0.5*30K", ylab="original")

LR.sum <- summary(lm(y ~ x))
abline(a = LR.sum$coefficients[1,1], b = LR.sum$coefficients[2,1], col="black")

# y <- bbb[, 2]; x <- bbb[, 1]
