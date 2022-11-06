rm(list=ls())
library(dplyr)
library(tidyr)

filepath <- "/projects/rpci/qzhu/users-qzhu/grant/UACA_R01/ELISA/"
PT.genotype <- read.table(paste0(filepath, "5th_ELISA.select_samples.txt"), header = TRUE, sep = "\t")

out.PT.genotype <- PT.genotype %>% select(SubjectID, rs6494889)

dat.450_570 <- read.table(paste0(filepath, "5thELISA_450nm-570nm"), header = TRUE, colClasses = c("numeric", "numeric", "numeric"))
dat.450_540 <- read.table(paste0(filepath, "5thELISA_450nm-540nm"), header = TRUE, colClasses = c("numeric", "numeric", "numeric"))


ID.info <- readxl::read_excel("/projects/rpci/qzhu/users-qzhu/grant/UACA_R01/ELISA/DistSamplesList_nowaka_Sep-16-2022-11-25-05.xls") # collection_alias
ID.info <- ID.info %>% select(Subjectid, `Collection Alias`) 
colnames(ID.info) <- c("Subjectid", "ID_alias")
ID.info <- unique( ID.info )
ID.info$ID_alias <- as.numeric(ID.info$ID_alias)
# patient "PT-00280540" has two records
# PT-00280540	1002201
# PT-00280540	130431113


dat.450_540 <- left_join(dat.450_540, ID.info, by=c("ID"="ID_alias"))  
dat.450_570 <- left_join(dat.450_570, ID.info, by=c("ID"="ID_alias"))
length(unique(dat.450_570$Subjectid)) # 39 PTs
length(unique(dat.450_570$ID)) # 40 records

dat.450_540.genotype <- unique(left_join(dat.450_540, out.PT.genotype, by=c("Subjectid"="SubjectID")))
dat.450_570.genotype <- unique(left_join(dat.450_570, out.PT.genotype, by=c("Subjectid"="SubjectID")))

# scatterplot

get_LR <- function(x, y, group, file="450-540nm") {
  print("lm coef:")
  print(summary(lm(y ~ x)))
  print("person:")
  print(cor(x, y, method = c("pearson")))
  print("spearman:")
  print(cor(x, y, method = c("spearman")))
  
  colors <- c("black", # Orange
              "deepskyblue2", # Light green
              "goldenrod2") # Darker green
  
  setwd("/user/xinyugao/ELISA_project/")
  pdf(paste0("5th_scatterplot_original.VS.30K_", file, ".pdf"))
  plot(x, y, pch=19,  xlim = c(min(x)-1, max(x)+1), ylim = c(min(y)-5, max(y)+5) ,col = colors[factor(group)],
       xlab="Concentration (original)", ylab="Concentration (30K)",
       main = paste0("scatterplot of original and 30K concentrations (", file, ")"))
  abline(lm(y ~ x), col = "red")
  legend("topright",
         legend = c("genotype=0", "genotype=1", "genotype=2"),
         # legend = levels(factor(group)),
         pch = 19,
         col = colors)
  dev.off()
}

get_LR(x = dat.450_540.genotype$original, y = dat.450_540.genotype$X30K, group=dat.450_540.genotype$rs6494889, file = "450-540nm")
get_LR(x = dat.450_570.genotype$original, y = dat.450_570.genotype$X30K, group=dat.450_570.genotype$rs6494889, file = "450-570nm")


get_LR_geno <- function(x, y, geno, file="450-540nm") {
  print("lm coef:")
  print(summary(lm(y ~ x)))
  print("person:")
  print(cor(x, y, method = c("pearson")))
  print("spearman:")
  print(cor(x, y, method = c("spearman")))

  
  setwd("/user/xinyugao/ELISA_project/")
  pdf(paste0("5th_scatterplot_original.VS.30K_genotype.", geno, "_",file, ".pdf"))
  plot(x, y, pch=19,  xlim = c(min(x)-1, max(x)+1), ylim = c(min(y)-5, max(y)+5) ,col = "black",
       xlab="Concentration (original)", ylab="Concentration (30K)",
       main = paste0("Original and 30K concentrations, genotype=", geno ," (", file, ")"))
  abline(lm(y ~ x), col = "red")
  dev.off()
}
# 540
dat.450_540.genotype.0 <- dat.450_540.genotype[dat.450_540.genotype$rs6494889 == 0,]
dat.450_540.genotype.1 <- dat.450_540.genotype[dat.450_540.genotype$rs6494889 == 1,]
# 570
dat.450_570.genotype.0 <- dat.450_570.genotype[dat.450_570.genotype$rs6494889 == 0,]
dat.450_570.genotype.1 <- dat.450_570.genotype[dat.450_570.genotype$rs6494889 == 1,]

get_LR_geno(x = dat.450_540.genotype.0$original, y = dat.450_540.genotype.0$X30K, geno = 0, file = "450-540nm"  )
get_LR_geno(x = dat.450_540.genotype.1$original, y = dat.450_540.genotype.1$X30K, geno = 1, file = "450-540nm"  )

get_LR_geno(x = dat.450_570.genotype.0$original, y = dat.450_570.genotype.0$X30K, geno = 0, file = "450-570nm"  )
get_LR_geno(x = dat.450_570.genotype.1$original, y = dat.450_570.genotype.1$X30K, geno = 1, file = "450-570nm"  )
#-----------------------------------------------------------------------------------------------------------------------------------------
# 450nm-570nm/450nm-540nm * (original/30K) (4 in total)
#  boxplot 3 groups (0,1,2) (option do not show outlier) + points (one point corresponds to a sample) + t-test (1 vs 0, 2 vs 0)


names <- c(rep("genotype=0", 19) , rep("genotype=1", 19) , rep("genotype=2", 2))
value.original <- c( dat.450_540.genotype$original[dat.450_540.genotype$rs6494889==0] , 
                     dat.450_540.genotype$original[dat.450_540.genotype$rs6494889==1],
                     dat.450_540.genotype$original[dat.450_540.genotype$rs6494889==2])
value.30K<- c( dat.450_540.genotype$X30K[dat.450_540.genotype$rs6494889==0] , 
                     dat.450_540.genotype$X30K[dat.450_540.genotype$rs6494889==1],
                     dat.450_540.genotype$X30K[dat.450_540.genotype$rs6494889==2])

data <- data.frame(names,value.original, value.30K)

##########################################################################################
# Basic boxplot
# 450-540 original
setwd("/user/xinyugao/ELISA_project/")
pdf(paste0("5th_boxplot_genotypes_450_540.original.pdf"))
boxplot(data$value.original ~ data$names , col=terrain.colors(4) , ylim=c(0, max(data$value.original) + 40),outline=FALSE,
        xlab = "Genotype Groups", ylab = "Concentration")

# Add data points
mylevels <- levels(factor(data$names))
levelProportions <- table(data$names)/nrow(data)
for(i in 1:length(mylevels)){
  
  thislevel <- mylevels[i]
  thisvalues <- data$value.original[data$names==thislevel]
  
  # take the x-axis indices and add a jitter, proportional to the N in each level
  myjitter <- jitter(rep(i, length(thisvalues)), amount=levelProportions[i]/2)
  points(myjitter, thisvalues, pch=20, col=rgb(0,0,0,.9)) 
}
legend("topright",
       legend = c("genotype=0 (n=19)", "genotype=1 (n=19)", "genotype=2 (n=2)"),
       fill = terrain.colors(4))
text(1, 190, "t-test for '1' and '0', p=0.69")
text(1, 200, "t-test for '2' and '0', p=0.78")
dev.off()
# 1 vs 0
t.test(data$value.original[data$names=="genotype=0"], data$value.original[data$names=="genotype=1"]) # p-value = 0.6883
# 2 vs 0
t.test(data$value.original[data$names=="genotype=0"], data$value.original[data$names=="genotype=2"]) # p-value = 0.7843

# 450-540 30K
setwd("/user/xinyugao/ELISA_project/")
pdf(paste0("5th_boxplot_genotypes_450_540.30K.pdf"))
boxplot(data$value.30K ~ data$names , col=terrain.colors(4) , ylim=c(0, max(data$value.30K) + 10),outline=FALSE,
        xlab = "Genotype Groups", ylab = "Concentration")

# Add data points
mylevels <- levels(factor(data$names))
levelProportions <- table(data$names)/nrow(data)
for(i in 1:length(mylevels)){
  
  thislevel <- mylevels[i]
  thisvalues <- data$value.30K[data$names==thislevel]
  
  # take the x-axis indices and add a jitter, proportional to the N in each level
  myjitter <- jitter(rep(i, length(thisvalues)), amount=levelProportions[i]/2)
  points(myjitter, thisvalues, pch=20, col=rgb(0,0,0,.9)) 
}
legend("topright",
       legend = c("genotype=0 (n=19)", "genotype=1 (n=19)", "genotype=2 (n=2)"),
       fill = terrain.colors(4))
text(1, 190, "t-test for '1' and '0', p=0.64")
text(1, 200, "t-test for '2' and '0', p=0.75")
dev.off()
# 1 vs 0
t.test(data$value.30K[data$names=="genotype=0"], data$value.30K[data$names=="genotype=1"]) # p-value = 0.6366
# 2 vs 0
t.test(data$value.30K[data$names=="genotype=0"], data$value.30K[data$names=="genotype=2"]) # p-value = 0.7538

#########################################################################################

# 450-570
names <- c(rep("genotype=0", 19) , rep("genotype=1", 19) , rep("genotype=2", 2))
value.original <- c( dat.450_570.genotype$original[dat.450_570.genotype$rs6494889==0] , 
                     dat.450_570.genotype$original[dat.450_570.genotype$rs6494889==1],
                     dat.450_570.genotype$original[dat.450_570.genotype$rs6494889==2])
value.30K<- c( dat.450_570.genotype$X30K[dat.450_570.genotype$rs6494889==0] , 
               dat.450_570.genotype$X30K[dat.450_570.genotype$rs6494889==1],
               dat.450_570.genotype$X30K[dat.450_570.genotype$rs6494889==2])

data <- data.frame(names,value.original, value.30K)

# 450-570 original
setwd("/user/xinyugao/ELISA_project/")
pdf(paste0("5th_boxplot_genotypes_450_570.original.pdf"))
boxplot(data$value.original ~ data$names , col=terrain.colors(4) , ylim=c(0, max(data$value.original) + 40),outline=FALSE,
        xlab = "Genotype Groups", ylab = "Concentration")

# Add data points
mylevels <- levels(factor(data$names))
levelProportions <- table(data$names)/nrow(data)
for(i in 1:length(mylevels)){
  
  thislevel <- mylevels[i]
  thisvalues <- data$value.original[data$names==thislevel]
  
  # take the x-axis indices and add a jitter, proportional to the N in each level
  myjitter <- jitter(rep(i, length(thisvalues)), amount=levelProportions[i]/2)
  points(myjitter, thisvalues, pch=20, col=rgb(0,0,0,.9)) 
}
legend("topright",
       legend = c("genotype=0 (n=19)", "genotype=1 (n=19)", "genotype=2 (n=2)"),
       fill = terrain.colors(4))
text(1, 190, "t-test for '1' and '0', p=0.70")
text(1, 200, "t-test for '2' and '0', p=0.79")
dev.off()
# 1 vs 0
t.test(data$value.original[data$names=="genotype=0"], data$value.original[data$names=="genotype=1"]) # p-value = 0.7003
# 2 vs 0
t.test(data$value.original[data$names=="genotype=0"], data$value.original[data$names=="genotype=2"]) # p-value = 0.7917



# 450-570 30K
setwd("/user/xinyugao/ELISA_project/")
pdf(paste0("5th_boxplot_genotypes_450_570.30K.pdf"))
boxplot(data$value.30K ~ data$names , col=terrain.colors(4) , ylim=c(0, max(data$value.30K) + 10),outline=FALSE,
        xlab = "Genotype Groups", ylab = "Concentration")

# Add data points
mylevels <- levels(factor(data$names))
levelProportions <- table(data$names)/nrow(data)
for(i in 1:length(mylevels)){
  
  thislevel <- mylevels[i]
  thisvalues <- data$value.30K[data$names==thislevel]
  
  # take the x-axis indices and add a jitter, proportional to the N in each level
  myjitter <- jitter(rep(i, length(thisvalues)), amount=levelProportions[i]/2)
  points(myjitter, thisvalues, pch=20, col=rgb(0,0,0,.9)) 
}
legend("topright",
       legend = c("genotype=0 (n=19)", "genotype=1 (n=19)", "genotype=2 (n=2)"),
       fill = terrain.colors(4))
text(1, 190, "t-test for '1' and '0', p=0.63")
text(1, 200, "t-test for '2' and '0', p=0.74")
dev.off()
# 1 vs 0
t.test(data$value.30K[data$names=="genotype=0"], data$value.30K[data$names=="genotype=1"]) # p-value = 0.6259
# 2 vs 0
t.test(data$value.30K[data$names=="genotype=0"], data$value.30K[data$names=="genotype=2"]) # p-value = 0.7436
