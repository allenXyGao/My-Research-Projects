#!/bin/sh
module use /projects/rpci/modules
module load plink/1.07


input=/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/cleaned_dataset/cleaned_rare_variants
output=/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/subtypes_Tcell_Bcell
sampleListFiles=/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/subtypes_Tcell_Bcell
declare -a pop
declare -a cohort_names
declare -a subtypes
cohort_names=(cohort1 cohort2)
subtypes=(Tcell Bcell)
pop=(ALLcasecontrol_cohort1 ALLcasecontrol_cohort2)


for i in ${!pop[@]}; do
	for j in ${!subtypes[@]}; do
		plink --bfile $input/CLEAN.rare.${pop[i]} \
	      		--keep $sampleListFiles/${subtypes[j]}_Meta/ALLEVC2.${cohort_names[i]}.${subtypes[j]}.sampleList.txt \
              		--make-bed \
              		--out $output/${subtypes[j]}_Meta/CLEAN.rare.${pop[i]} \
	      		--noweb	
	done
done


