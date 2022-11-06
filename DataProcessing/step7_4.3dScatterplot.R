## Run this file directly in Rstudio

#######################################################
### 1) generate the removed samples list for each pop
######################################################
pop <- c("ALLcasecontrol_cohort1", "ALLcasecontrol_cohort2",
         "AMLMDScasecontrol_cohort1", "AMLMDScasecontrol_cohort2")

indicator_remove <- rep(FALSE, length(pop))
# outlier ind: eth file
for (i in 1:length(pop)){
  filepath <- paste("/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/EigenStrat_Results/",
                    pop[i], "/LDprune_autosome-pca/", sep="")
  pca1 <- read.table(file=paste(filepath, "LDprune_autosome-pca1.fam", sep=""))
  pca2 <- read.table(file=paste(filepath, "LDprune_autosome-pca2.fam", sep=""))
  remInd <- setdiff(pca1[,2], pca2[,2])
  if (length(remInd)>0){
    indicator_remove[i] <- TRUE
    remInd1 <- pca1[match(remInd, pca1[,2]), 1:2]
    write.table(remInd1, file=paste(filepath, "removedSamples.", pop[i], sep=""), sep="\t",
                row.names=FALSE, col.names=FALSE, quote=FALSE)
  }
}

pop <- c("ALLcasecontrol_cohort1", "ALLcasecontrol_cohort2",
         "AMLMDScasecontrol_cohort1", "AMLMDScasecontrol_cohort2")
for (i in 1:length(pop)){
  if (indicator_remove[i] == FALSE) {
    next;
  }
  s <- read.table(file=paste("/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/EigenStrat_Results/", pop[i],
                           "/LDprune_autosome-pca/removedSamples.", pop[i], sep=""), header=FALSE)
  s1 <- paste(s[,1], s[,2], sep=":")
  assign(paste("remove.", pop[i], sep=""), s1)
}

##################################################
### 2) calculate the R2 for eigenvalues
##################################################
pop <- c("ALLcasecontrol_cohort1", "ALLcasecontrol_cohort2",
         "AMLMDScasecontrol_cohort1", "AMLMDScasecontrol_cohort2")

r2 <- c()
for (i in 1:length(pop)){
  filepath <- paste("/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/EigenStrat_Results/", pop[i],
                    "/LDprune_autosome-pca/", sep="")
  twout <- read.table(file=paste(filepath, "LDprune_autosome-pca2.twout", sep=""))
  sig.n <- sum(twout[1:20,5]<=0.05)  ## use top 20 to calculate, sig.n<20
  print(sig.n)
  eigen.v <- read.table(file=paste(filepath, "LDprune_autosome-pca2.eval", sep=""))
  r2 <- c(r2, sum(eigen.v[1:sig.n, 1])/sum(eigen.v[,1]))
}
## > r2
## [1] 0.009503656 0.005723215 0.012417999 0.007224144

##################################################
### 3) check the PCA results using 3D scatterplot:
##################################################

pop <- c("ALLcasecontrol_cohort1", "ALLcasecontrol_cohort2",
         "AMLMDScasecontrol_cohort1", "AMLMDScasecontrol_cohort2")
library(rgl)
# for (i in 1:length(pop)){
#   X11(width=10,height=10)
#   dat1 <- read.table(file=paste("/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/EigenStrat_Results/",
#                                 pop[i], "/LDprune_autosome-pca1.eth.out", sep=""), header=T, as.is=T)
#   dat2 <- read.table(file=paste("/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/EigenStrat_Results/",
#                                 pop[i], "/LDprune_autosome-pca2.eth.out", sep=""), header=T, as.is=T)
#   xrange <- round(range(dat1[,3]), 1)
#   yrange <- round(range(dat1[,4]), 1)
#   zrange <- round(range(dat1[,5]), 1)
#   ## plot using the first 3 eigenvectors
#   ## 1) plot pca1 by marking those removed in pca2
#   clr <- rep("black", length=nrow(dat1))
#   if (indicator_remove[i] == TRUE) {
#     clr[match(eval(parse(text=paste("remove.", pop[i], sep=""))), dat1[,1])] <- "red"
#   }
#   plot3d(dat1[,3], dat1[,4], dat1[,5], size=3, col=clr, xlab="PCA1", ylab="PCA2", zlab="PCA3", xlim=xrange, ylim=yrange, zlim=zrange)
#   rgl.postscript(file=paste("/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/EigenStrat_Results/",
#                             pop[i], "/LDprune_autosome-pca/", pop[i], ".pca1.pdf", sep=""), fmt="pdf")
#   ## 2) plot pca2, using same x,y,z limit.
#   plot3d(dat2[,3], dat2[,4], dat2[,5], size=3, xlab="PCA1", ylab="PCA2", zlab="PCA3", xlim=xrange, ylim=yrange, zlim=zrange)
#   rgl.postscript(file=paste("/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/EigenStrat_Results/",
#                             pop[i], "/LDprune_autosome-pca/", pop[i], ".pca2.pdf", sep=""), fmt="pdf")
# }


##########################################################################################
this_3dplot.pc1 <- function(i) {
  dat1 <- read.table(file=paste("/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/EigenStrat_Results/",
                                pop[i], "/LDprune_autosome-pca1.eth.out", sep=""), header=T, as.is=T)
  dat2 <- read.table(file=paste("/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/EigenStrat_Results/",
                                pop[i], "/LDprune_autosome-pca2.eth.out", sep=""), header=T, as.is=T)
  xrange <- round(range(dat1[,3]), 1)
  yrange <- round(range(dat1[,4]), 1)
  zrange <- round(range(dat1[,5]), 1)
  ## plot using the first 3 eigenvectors
  ## 1) plot pca1 by marking those removed in pca2
  clr <- rep("black", length=nrow(dat1))
  if (indicator_remove[i] == TRUE) {
    clr[match(eval(parse(text=paste("remove.", pop[i], sep=""))), dat1[,1])] <- "red"
  }
  plot3d(dat1[,3], dat1[,4], dat1[,5], size=3, col=clr, xlab="PCA1", ylab="PCA2", zlab="PCA3", xlim=xrange, ylim=yrange, zlim=zrange)
  # bgplot3d({
  #   plot.new()
  #   title(main=paste(pop[i], "pc1 data"), line=3)
  # })
}

this_3dplot.pc2 <- function(i) {
  dat1 <- read.table(file=paste("/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/EigenStrat_Results/",
                                pop[i], "/LDprune_autosome-pca1.eth.out", sep=""), header=T, as.is=T)
  dat2 <- read.table(file=paste("/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/EigenStrat_Results/",
                                pop[i], "/LDprune_autosome-pca2.eth.out", sep=""), header=T, as.is=T)
  xrange <- round(range(dat2[,3])*2, 1)
  yrange <- round(range(dat2[,4])*2, 1)
  zrange <- round(range(dat2[,5])*2, 1)
  plot3d(dat2[,3], dat2[,4], dat2[,5], size=3, col="black",xlab="PCA1", ylab="PCA2", zlab="PCA3", xlim=xrange, ylim=yrange, zlim=zrange)
  # bgplot3d({
  #   plot.new()
  #   title(main=paste(pop[i], "pc2 data"), line=3)
  # })
}


# ALLcasecontrol_cohort1
this_3dplot.pc1(1)
this_3dplot.pc2(1)

# ALLcasecontrol_cohort2
this_3dplot.pc1(2)
this_3dplot.pc2(2)

# AMLMDScasecontrol_cohort1
this_3dplot.pc1(3)
this_3dplot.pc2(3)

# AMLMDScasecontrol_cohort2
this_3dplot.pc1(4)
this_3dplot.pc2(4)


