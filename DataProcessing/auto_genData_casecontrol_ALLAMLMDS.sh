#!/bin/sh

input=/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/input
output=/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output

module use /projects/rpci/modules
module load plink/1.07

caseFi=/projects/rpci/qzhu/eschille/transplant.GRANT/data/Recipients
ctrlFi=/projects/rpci/qzhu/qzhu/mywork/LaraBMTexome/input

declare -a cohorts
declare -a casetypes
declare -a cohort_names
declare -a casetype_names
declare -a control_cohorts
cohorts=(EAcohort1 EAcohort2)
cohort_names=(cohort1 cohort2)
casetypes=(ALL.samples.txt AML.MDS.samples.txt)
casetype_names=(ALLcasecontrol AMLMDScasecontrol)
control_cohorts=(BMT.snpQC.final.donor.EA.cohort1.fam BMT.snpQC.final.donor.EA.cohort2.fam)

for i in ${!cohorts[@]}; do
	for j in ${!casetypes[@]}; do
		awk '{print $1"\t"$2}' $caseFi/${cohorts[i]}/${casetypes[j]} > $input/samplelist.${cohort_names[i]}.${casetype_names[j]}
		awk '{print $1"\t"$2}' $ctrlFi/${control_cohorts[i]} >> $input/samplelist.${cohort_names[i]}.${casetype_names[j]}
		plink --bfile /projects/rpci/qzhu/qliu7/mywork/20140609.BMT_exome/output/otherCheck/merge_cleanedAllSamples \
			--keep $input/samplelist.${cohort_names[i]}.${casetype_names[j]} \
			--geno 0.02 --hwe 1e-6 --mind 0.02 \
			--make-bed \
			--out $output/${casetype_names[j]}_${cohort_names[i]} \
			--noweb
	done
done





#### cases
# ALL (cases in cohort 1): /projects/rpci/qzhu/eschille/transplant.GRANT/data/Recipients/EAcohort1/ALL.samples.txt
# AML.MDS(cases in cohort 1): /projects/rpci/qzhu/eschille/transplant.GRANT/data/Recipients/EAcohort1/AML.MDS.samples.txt
# ALL(cases in cohort 2): /projects/rpci/qzhu/eschille/transplant.GRANT/data/Recipients/EAcohort2/ALL.samples.txt
# AML.MDS(cases in cohort 2): /projects/rpci/qzhu/eschille/transplant.GRANT/data/Recipients/EAcohort2/AML.MDS.samples.txt

#### control
# (control in cohort 1): /projects/rpci/qzhu/qzhu/mywork/LaraBMTexome/input/BMT.snpQC.final.donor.EA.cohort1.fam
# (control in cohort 2): /projects/rpci/qzhu/qzhu/mywork/LaraBMTexome/input/BMT.snpQC.final.donor.EA.cohort2.fam





