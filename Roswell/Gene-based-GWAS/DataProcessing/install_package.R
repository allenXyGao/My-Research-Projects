# install skatMeta and survivalSKAT and their dependencies
# install.packages("coxme", lib="/projects/rpci/qzhu/xinyu/Software/R/R-4.0.2", dependencies = TRUE, repos="http://lib.stat.cmu.edu/R/CRAN")
# install.packages("CompQuadForm", lib="/projects/rpci/qzhu/xinyu/Software/R/R-4.0.2", dependencies = TRUE, repos="http://lib.stat.cmu.edu/R/CRAN")
# install.packages("cmprsk", lib="/projects/rpci/qzhu/xinyu/Software/R/R-4.0.2", dependencies = TRUE, repos="http://lib.stat.cmu.edu/R/CRAN")
# install.packages("/projects/rpci/qzhu/qzhu/software/R/skatMeta_1.4.3.tar.gz", repos = NULL, type="source", 
#                  lib="/projects/rpci/qzhu/xinyu/Software/R/R-4.0.2")
# install.packages("/projects/rpci/qzhu/qzhu/software/R/survivalSKAT_0.3.0000.tar.gz", repos = NULL, type="source", 
#                  lib="/projects/rpci/qzhu/xinyu/Software/R/R-4.0.2")

#library(coxme, lib.loc = "/projects/rpci/qzhu/xinyu/Software/R/R-4.0.2")


# install.packages("package_name", lib="/projects/rpci/qzhu/xinyu/Software/R/R-4.0.2")   # lib: the path where you store R packages
mypaths <- .libPaths()
mypaths <- c(mypaths, "/projects/rpci/qzhu/xinyu/Software/R/R-4.0.2")
.libPaths(mypaths)

library(skatMeta)
library(survivalSKAT)

head(snpInfo.functional)
# output -> txt 
# 1st column: Genes
# 2nd column: Name
temp <- split(snpInfo.functional, f=snpInfo.functional$Genes)
head(temp)

temp.data <- snpInfo.functional[, colnames(snpInfo.functional)[c(5,1)]]
head(temp.data)
filepath="/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/SSD_results/"
write.table(temp.data, file = paste0(filepath, "File.SetID"), 
            row.names = FALSE, col.names = FALSE, 
            quote = FALSE, sep=" ")


library(SKAT)
####################################################################################################################



####################################################################################################################
# install package seqMeta
# install.packages("/projects/rpci/qzhu/xinyu/Software/seqMeta-1.6.7.tar.gz", repos = NULL, type="source", 
#                   lib="/projects/rpci/qzhu/xinyu/Software/R/R-4.0.2")
mypaths <- .libPaths()
mypaths <- c(mypaths, "/projects/rpci/qzhu/xinyu/Software/R/R-4.0.2")
.libPaths(mypaths)

library(seqMeta)


####################################################################################################################
# install package MetaSKAT
# install.packages("/projects/rpci/qzhu/xinyu/Software/MetaSKAT_0.81.tar.gz", repos = NULL, type="source", 
#                    lib="/projects/rpci/qzhu/xinyu/Software/R/R-4.0.2")
library(MetaSKAT)











