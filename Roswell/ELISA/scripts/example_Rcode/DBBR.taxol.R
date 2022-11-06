# based on "11-19.3" in /projects/rpci/qzhu/eschille/BreastGWAS/No_Y_Chrom/maf1.noHM/DBBR/Scripts/dbbr.associations.R
library(survival)

model = function(cov, geno) {
	res = coxph(Surv(days_os, PatientStatus_Description) ~ DxAge + BMI2 + as.factor(ER) + as.factor(PR) + as.factor(HER2) + as.factor(Hormonal.Therapy) +
                as.factor(Radiation.Therapy) + as.factor(Grade_Description) + as.factor(StageNEW) + as.factor(Surgery.Type) + geno, cbind(cov, geno=geno)) 
	return(summary(res)$coefficients["geno",])	
}

# no HER2 as covariate
model2 = function(cov, geno) {
        res = coxph(Surv(days_os, PatientStatus_Description) ~ DxAge + BMI2 + as.factor(ER) + as.factor(PR) + as.factor(Hormonal.Therapy) +
                as.factor(Radiation.Therapy) + as.factor(Grade_Description) + as.factor(StageNEW) + as.factor(Surgery.Type) + geno, cbind(cov, geno=geno))
        return(summary(res)$coefficients["geno",])
}

pheno = read.csv("/projects/rpci/qzhu/eschille/BreastGWAS/No_Y_Chrom/maf1.noHM/DBBR/Input/DBBR.pheno.csv", na.strings = c("NA", ''))
pheno = within(pheno, ER <- relevel(as.factor(ER), ref = "2"))
pheno = within(pheno, PR <- relevel(as.factor(PR), ref = "2"))
pheno = within(pheno, HER2 <- relevel(as.factor(HER2), ref = "2"))
pheno$DxAge = as.numeric(as.character(pheno$DxAge))
#pheno = within(pheno, Surgery.Type.clean <- relevel(as.factor(Surgery.Type.clean), ref = "None"))
pheno = within(pheno, Surgery.Type <- relevel(as.factor(Surgery.Type), ref = "None"))

pheno$rs11855431 = as.numeric(as.character(pheno$rs11855431))
pheno$rs6494889 = as.numeric(as.character(pheno$rs6494889))
pheno$rs28607477 = as.numeric(as.character(pheno$rs28607477))
pheno$rs720251 = as.numeric(as.character(pheno$rs720251))

#doxo.yes = pheno$SubjectID[!is.na(pheno$doxo.before.surgery) & pheno$doxo.before.surgery=="YES"]
#anti.yes = pheno$SubjectID[!is.na(pheno$anti.before.surgery) & pheno$anti.before.surgery=="YES"]
#yeses = c(doxo.yes, anti.yes)
#yeses = unique(yeses)

taxol.pts = read.table("/projects/rpci/qzhu/qzhu/mywork/PathwaysGWAS/input/DBBR.taxol.pts", header=T, sep="\t")
ind = match(pheno$SubjectID, taxol.pts$SubjectID)
pheno$taxol = !is.na(ind)

# pts received anthracycline or HER2 target therapy
newdat = pheno[(!is.na(pheno$Anthracyline.clean) & pheno$Anthracyline.clean=="YES") | (!is.na(pheno$anti.Her2.clean ) & !is.na(pheno$HER2) & pheno$HER2=="1" & pheno$anti.Her2.clean=="YES"),]
print(nrow(newdat))

cov = newdat[,c("days_os", "PatientStatus_Description", "DxAge", "BMI2", "ER", "PR", "HER2", "Hormonal.Therapy", "Radiation.Therapy", "Grade_Description", "StageNEW", "Surgery.Type")]
geno = newdat[,c("rs6494889", "rs28607477", "rs720251", "rs11855431")]
res = apply(geno, 2, function(x){model(cov=cov, geno=x)})
write.table(t(res), file="/projects/rpci/qzhu/qzhu/mywork/PathwaysGWAS/output/DBBR/DBBR.treatmentsubgroup.test", append=F, quote=F, sep="\t", col.names=T, row.names=T)	#n= 396, number of events= 63

# pts without taxol
tmp = newdat[!newdat$taxol,]
print(nrow(tmp))
tmp$Grade_Description[tmp$Grade_Description==1] = 2 # only 7 pts have grade "1"; got warning (Loglik converged before variable  8,9 ; coefficient may be infinite) when fitting model without merging grade 1 and 2
cov = tmp[,c("days_os", "PatientStatus_Description", "DxAge", "BMI2", "ER", "PR", "HER2", "Hormonal.Therapy", "Radiation.Therapy", "Grade_Description", "StageNEW", "Surgery.Type")]
geno = tmp[,c("rs6494889", "rs28607477", "rs720251", "rs11855431")]
res = apply(geno, 2, function(x){model(cov=cov, geno=x)})
write.table(t(res), file="/projects/rpci/qzhu/qzhu/mywork/PathwaysGWAS/output/DBBR/DBBR.treatmentsubgroupWOtaxol.test", append=F, quote=F, sep="\t", col.names=T, row.names=T)	# n= 101, number of events= 16 

# pts with taxol
tmp = newdat[newdat$taxol,]
print(nrow(tmp))
cov = tmp[,c("days_os", "PatientStatus_Description", "DxAge", "BMI2", "ER", "PR", "HER2", "Hormonal.Therapy", "Radiation.Therapy", "Grade_Description", "StageNEW", "Surgery.Type")]
geno = tmp[,c("rs6494889", "rs28607477", "rs720251", "rs11855431")]
res = apply(geno, 2, function(x){model(cov=cov, geno=x)})
write.table(t(res), file="/projects/rpci/qzhu/qzhu/mywork/PathwaysGWAS/output/DBBR/DBBR.treatmentsubgroupWtaxol.test", append=F, quote=F, sep="\t", col.names=T, row.names=T)	# n= 295, number of events= 47

# pts received anthracycline
newdat = pheno[(!is.na(pheno$Anthracyline.clean) & pheno$Anthracyline.clean=="YES"),]
print(nrow(newdat))
cov = newdat[,c("days_os", "PatientStatus_Description", "DxAge", "BMI2", "ER", "PR", "HER2", "Hormonal.Therapy", "Radiation.Therapy", "Grade_Description", "StageNEW", "Surgery.Type")]
geno = newdat[,c("rs6494889", "rs28607477", "rs720251", "rs11855431")]
res = apply(geno, 2, function(x){model(cov=cov, geno=x)})
write.table(t(res), file="/projects/rpci/qzhu/qzhu/mywork/PathwaysGWAS/output/DBBR/DBBR.Anthracyline.test", append=F, quote=F, sep="\t", col.names=T, row.names=T)

# pts received anthracycline but without taxol
tmp = newdat[!newdat$taxol,]
print(nrow(tmp))
tmp$Grade_Description[tmp$Grade_Description==1] = 2 # only 5 pts have grade "1"; got warning (Loglik converged before variable  8,9 ; coefficient may be infinite) when fitting model without merging grade 1 and 2
#tmp$StageNEW[tmp$StageNEW=="I"] = "II" # 17 pts have stage "1"; got warning (Loglik converged before variable  8,9 ; coefficient may be infinite) when fitting model without merging
tmp$StageNEW[tmp$StageNEW=="III/IV"] = "II" # 17 pts have stage "1"; got warning (Loglik converged before variable  8,9 ; coefficient may be infinite) when fitting model without merging
cov = tmp[,c("days_os", "PatientStatus_Description", "DxAge", "BMI2", "ER", "PR", "HER2", "Hormonal.Therapy", "Radiation.Therapy", "Grade_Description", "StageNEW", "Surgery.Type")]
geno = tmp[,c("rs6494889", "rs28607477", "rs720251", "rs11855431")]
res = apply(geno, 2, function(x){model2(cov=cov, geno=x)})	# have to remove HER2 from covariates as only 5 pts are HER2+
write.table(t(res), file="/projects/rpci/qzhu/qzhu/mywork/PathwaysGWAS/output/DBBR/DBBR.AnthracylineWOtaxol.test", append=F, quote=F, sep="\t", col.names=T, row.names=T)

# pts received anthracycline and taxol
tmp = newdat[newdat$taxol,]
print(nrow(tmp))
cov = tmp[,c("days_os", "PatientStatus_Description", "DxAge", "BMI2", "ER", "PR", "HER2", "Hormonal.Therapy", "Radiation.Therapy", "Grade_Description", "StageNEW", "Surgery.Type")]
geno = tmp[,c("rs6494889", "rs28607477", "rs720251", "rs11855431")]
res = apply(geno, 2, function(x){model(cov=cov, geno=x)})
write.table(t(res), file="/projects/rpci/qzhu/qzhu/mywork/PathwaysGWAS/output/DBBR/DBBR.AnthracylineWtaxol.test", append=F, quote=F, sep="\t", col.names=T, row.names=T)



