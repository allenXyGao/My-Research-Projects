library(SKAT)

# load simplified status for DNMT3A
pop <- "AMLMDScasecontrol_cohort2"
gene.name <- "DNMT3A"
filepath <- "/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/variant_plots/"
info.status.cohort1 <- read.table(file = paste0(filepath, gene.name, "_info.status_cohort1.txt"), 
           header = TRUE, colClasses = c(rep("character", 2)))
info.status.cohort2 <- read.table(file = paste0(filepath, gene.name, "_info.status_cohort2.txt"), 
                                  header = TRUE, colClasses = c(rep("character", 2)))
info.status.cohort1
info.status.cohort2



getDiseaseCohortNames <- function(dataFile="ALLcasecontrol_cohort1") {
  split.Full_name <- unlist(strsplit(dataFile, split="_"))
  cohort.name <- split.Full_name[2]
  disease.name <- "AMLMDS"
  if (paste(unlist(strsplit(split.Full_name[1], split =""))[1:3], collapse = "") == "ALL") {
    disease.name <- "ALL"
  }
  return(list("disease_name"=disease.name, "cohort_name"=cohort.name))
}


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

#case.control_status[!is.na(match(case.control_status$ID, info.status.cohort1$ID)),]


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

change_col <- function(col.mat, info.status) {
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
        # further check if status is remission or relapse or NA
        ID <- case.control_status$ID[inds.ordered[matched_sample]]
        status <- info.status$status[info.status$ID == ID]
        if (is.na(status)) {
          col.mat[i, matched_sample] = "red"
        }
        else if (status == "CR") {
          col.mat[i, matched_sample] = "orange"
        }
        else if (status == "Relapse") {
          col.mat[i, matched_sample] = "hotpink"
        }
      }
      else {
        col.mat[i, matched_sample] = "blue"
      }
    }
  }
  return(col.mat)
}

#col.mat <- change_col(col.mat = col.mat, info.status = info.status.cohort1)
col.mat <- change_col(col.mat = col.mat, info.status = info.status.cohort2)
col.mat <- as.matrix(col.mat, nrow = n_snps, ncol = ninds)

## need reorder column and rows
##### reorder rows
if (dim(col.mat)[1] == 1) {
  col.mat <- col.mat
} else {
  row.sums.raw <- rowSums(col.mat != "grey")
  row.order.info <- sort(row.sums.raw, decreasing = TRUE,index.return=TRUE)
  col.mat <- col.mat[row.order.info$ix, ] # update col.mat
  snps <- snps[row.order.info$ix]# update snp names
}


#### reorder columns
# each column can have only one color -> identify the category of the first non-grey element  
carriers_ID <- case.control_status$ID[inds.ordered]
column.status.info <- rep("grey", ninds)
for (j in 1:ninds) {
  match.info <- col.mat[,j][col.mat[, j] != "grey"]
  if (is.na(match.info)) {next;}
  column.status.info[j] <- match.info[1]
}
#column.status.info
column.status.info <- factor(column.status.info, levels = c("orange", "hotpink", "red", "blue"))
column.order.index <- order(column.status.info)
column.status.info[column.order.index]

col.mat <- col.mat[, column.order.index] # update col.mat
carriers_ID <- carriers_ID[column.order.index] # update carriers ID
col.mat <- matrix(col.mat, nrow = n_snps, ncol = ninds)

###############################
# step 4: plot
setwd("/user/xinyugao/SKATO_Project/plots_carriers_variants/")
pdf(paste0("CR_NonCR_summary_", pop, "_",gene.name, ".pdf"))
xleft = 0:(ninds-1)
xright = xleft+1
vertical.len <- 0.25
plot(0, 0, type = "n", xlim = c(0, ninds), ylim = c(0, n_snps+1), 
     axes = FALSE, xlab = "", ylab = "" )
rect(xleft, n_snps/(1/vertical.len), xright, n_snps/(1/vertical.len)+vertical.len, border = "light gray" , col=case.control.col)
axis(2, at=n_snps/(1/vertical.len)+vertical.len/2, labels="case/control status", tick=F, lty=6, 
     pos=0.3, las=2, cex.axis=0.6)
axis(3, at=xleft+0.5, labels=carriers_ID, tick=F , pos=n_snps/(1/vertical.len)+vertical.len, 
     las=2, cex.axis=0.5, xpd = TRUE)
start = n_snps / (1/vertical.len)
ybottom <- start - seq(vertical.len, n_snps / (1/vertical.len), by = vertical.len)
ytop <- ybottom + vertical.len

for (j in 1:n_snps) {
  rect(xleft, ybottom[j], xright, ytop[j], border="light grey", col=col.mat[j, ])
  axis(2, at=ybottom[j]+vertical.len/2, labels = snps[j], tick=F, lty=6, pos=0, las=2, cex.axis=0.6)
  axis(4, at=ybottom[j]+vertical.len/2, labels=sum(col.mat[j, ]!="grey"), 
       tick=F, lty=6, pos=ninds-0.2, las=2, cex.axis=0.7)
}
axis(1, at=xleft+0.5, labels=colSums(col.mat != "grey"), tick=F, pos=0, las=2, cex.axis=0.5, xpd = TRUE)
# legend(x=0, y=-0.2, legend=c("CR", "non CR", "case (NA)", "control"), 
#        pt.bg=c("orange", "hotpink", "red", "blue"), pch=22, cex=0.6, xpd = TRUE)
legend(x=0, y=-0.2, legend=c("CR", "non CR"), 
        pt.bg=c("orange", "hotpink"), pch=22, cex=0.6, xpd = TRUE)
dev.off()


######################################################################################