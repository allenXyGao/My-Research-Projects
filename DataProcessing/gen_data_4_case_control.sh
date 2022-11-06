#!/bin/sh

module use /projects/rpci/modules
module load plink/1.07

ctrlFi=/projects/rpci/qzhu/qzhu/mywork/LaraBMTexome/input/BMT.snpQC.final.donor.EA.cohort1.fam
caseFi=/projects/rpci/qzhu/eschille/transplant.GRANT/data/Recipients/EAcohort1/ALL.samples.txt

input=/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/input
output=/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output

#### cases
# ALL (cases in cohort 1): /projects/rpci/qzhu/eschille/transplant.GRANT/data/Recipients/EAcohort1/ALL.samples.txt

# AML.MDS(cases in cohort 1): /projects/rpci/qzhu/eschille/transplant.GRANT/data/Recipients/EAcohort1/AML.MDS.samples.txt

# ALL(cases in cohort 2): /projects/rpci/qzhu/eschille/transplant.GRANT/data/Recipients/EAcohort2/ALL.samples.txt

# AML.MDS(cases in cohort 2): /projects/rpci/qzhu/eschille/transplant.GRANT/data/Recipients/EAcohort2/AML.MDS.samples.txt


#### control

# (control in cohort 1): /projects/rpci/qzhu/qzhu/mywork/LaraBMTexome/input/BMT.snpQC.final.donor.EA.cohort1.fam

# (control in cohort 2): /projects/rpci/qzhu/qzhu/mywork/LaraBMTexome/input/BMT.snpQC.final.donor.EA.cohort2.fam



awk '{print $1"\t"$2}' $caseFi > $input/samplelist
awk '{print $1"\t"$2}' $ctrlFi >> $input/samplelist


plink --bfile /projects/rpci/qzhu/qliu7/mywork/20140609.BMT_exome/output/otherCheck/merge_cleanedAllSamples --keep $input/samplelist \
	--geno 0.02 --hwe 1e-6 --mind 0.02 \
	--make-bed \
	--out $output/ALLcasecontrol.cohort1 \
	--noweb



