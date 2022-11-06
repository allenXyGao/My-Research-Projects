### plot case/control carriers of variants


library(SKAT)
getDiseaseCohortNames <- function(dataFile="ALLcasecontrol_cohort1") {
  split.Full_name <- unlist(strsplit(dataFile, split="_"))
  cohort.name <- split.Full_name[2]
  disease.name <- "AMLMDS"
  if (paste(unlist(strsplit(split.Full_name[1], split =""))[1:3], collapse = "") == "ALL") {
    disease.name <- "ALL"
  }
  return(list("disease_name"=disease.name, "cohort_name"=cohort.name))
}

# pop <- "ALLcasecontrol_cohort2"
# gene.name <- "EVC2"
pop <- "AMLMDScasecontrol_cohort2"
gene.name <- "DNMT3A"

names_case_cohort <- getDiseaseCohortNames(dataFile = pop)
filepath="/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/SSD_results/"
res.open_ssd <- Open_SSD(File.SSD = paste0(filepath, names_case_cohort$disease_name,"casecontrol_", names_case_cohort$cohort_name,"/", "File.SSD"),
                         File.Info = paste0(filepath, names_case_cohort$disease_name,"casecontrol_", names_case_cohort$cohort_name,"/","File.Info"))

filepath.model <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/NullModels_OrderedPheno/"
load(file=paste0(filepath.model, "nullModel_orderedPheno_", pop, ".RData"))
phenotype_data <- res.phenotype_nullModel$phenotype_Data
pheno_with_FIDIID <- res.phenotype_nullModel$pheno_with_FIDIID


######################################
# step 1: list variants
index <- res.open_ssd$SetInfo[res.open_ssd$SetInfo$SetID == gene.name, "SetIndex"]
geno.mat <- Get_Genotypes_SSD(SSD_INFO = res.open_ssd,
                              Set_Index = index, is_ID = TRUE)
snps <- colnames(geno.mat)

#####################################
# step 2: order samples by case/control status
case.control_status <- data.frame("ID"=pheno_with_FIDIID$IID, "status"=pheno_with_FIDIID$phenotype)
out <- strsplit(case.control_status$ID, "([A-Z]-)", perl = TRUE)
case.control_status$ID <- gsub("[-_]", "", do.call(rbind, out)[,2])
head(case.control_status) # 1: case; 0: control




###############################
# step 3: change colors
carriers <- apply(geno.mat, 2, function(x) which(x!=0 & x!=9, arr.ind = T))
inds <- unique(unlist(carriers))
inds.ordered <- inds[order(inds)]
# case/control
ind.case <- match(inds.ordered, which(case.control_status$status == 1))
ind.case = inds.ordered[!is.na(ind.case)]
ind.control <- match(inds.ordered, which(case.control_status$status == 0))
ind.control = inds.ordered[!is.na(ind.control)]
case.control.col <- c(rep("red", length(ind.case)), rep("blue", length(ind.control)))

n_snps <- length(snps); ninds <- length(inds)
col.mat <- matrix("grey", nrow = n_snps, ncol = ninds)

for (i in 1:n_snps) {
  print(paste0("current snp is: ", snps[i]))
  #
  if (is.na(carriers[snps[i]])) {
    cur_carriers <- carriers
  } else {
    cur_carriers <- carriers[[i]]
  }

  for (j in 1:length(cur_carriers)) {
    matched_sample <- which(!is.na(match(inds.ordered, cur_carriers[j])))
    if (matched_sample <= length(ind.case)) {
      col.mat[i, matched_sample] = "red"
    }
    else {
      col.mat[i, matched_sample] = "blue"
    }
  }
}

carriers_ID <- case.control_status$ID[inds.ordered]
res.variant_sample <- list()
res.variant_sample$carriers_ID <- carriers_ID
res.variant_sample$col.mat <- col.mat
res.variant_sample$inds.ordered <- inds.ordered
res.variant_sample$ind.case <- ind.case
res.variant_sample$ind.control <- ind.control
res.variant_sample$snps <- snps
filepath.out <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/variant_plots/"
save(res.variant_sample, file=paste0(filepath.out, pop, "_", gene.name,  "_colorMat.RData"))


###############################
# step 4: plot
setwd("/user/xinyugao/SKATO_Project/plots_carriers_variants/")
pdf(paste0("summary_", pop, "_",gene.name, ".pdf"))
carriers_ID <- case.control_status$ID[inds.ordered]
xleft = 0:(ninds-1)
xright = xleft+1
vertical.len <- 0.25
plot(0, 0, type = "n", xlim = c(0, ninds), ylim = c(0, n_snps+1), 
     axes = FALSE, xlab = "", ylab = "" )
rect(xleft, n_snps/(1/vertical.len), xright, n_snps/(1/vertical.len)+vertical.len, border = "light gray" , col=case.control.col)
axis(2, at=n_snps/(1/vertical.len)+vertical.len/2, labels="case/control status", tick=F, lty=6, 
     pos=0.3, las=2, cex.axis=0.6)
axis(3, at=xleft+0.5, labels=carriers_ID, tick=F , pos=n_snps/(1/vertical.len)+vertical.len, 
     las=2, cex.axis=0.7, xpd = TRUE)
start = n_snps / (1/vertical.len)
ybottom <- start - seq(vertical.len, n_snps / (1/vertical.len), by = vertical.len)
ytop <- ybottom + vertical.len

for (j in 1:n_snps) {
  rect(xleft, ybottom[j], xright, ytop[j], border="light grey", col=col.mat[j, ])
  axis(2, at=ybottom[j]+vertical.len/2, labels = snps[j], tick=F, lty=6, pos=0, las=2, cex.axis=0.6)
  axis(4, at=ybottom[j]+vertical.len/2, labels=sum(col.mat[j, ]!="grey"), 
       tick=F, lty=6, pos=ninds-0.2, las=2, cex.axis=0.7)
}
axis(1, at=xleft+0.5, labels=colSums(col.mat != "grey"), tick=F, pos=0, las=2, cex.axis=0.7, xpd = TRUE)
legend(x=0, y=-0.2, legend=c("case", "control"), 
       pt.bg=c("red", "blue", "grey"), pch=22, cex=0.6, xpd = TRUE)
dev.off()


######################################################################################