#!/bin/sh
input=/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/cleaned_dataset/cleaned_rare_maf005 
samplelist=/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/input
output=/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/cleaned_dataset


module use /projects/rpci/modules
module load plink/1.07


declare -a subtypes
declare -a pops
subtypes=(Tcell Bcell)
pops=(ALLcasecontrol_cohort1 ALLcasecontrol_cohort2)



for i in ${!pops[@]}; do
	for j in ${!subtypes[@]}; do

	        plink --bfile $input/CLEAN.rare005.${pops[i]} \
	      	      --keep $samplelist/samplelist.${subtypes[j]}/${pops[i]}/${pops[i]}.samplelist.${subtypes[j]}.txt \
              	      --make-bed \
              	      --out $output/cleaned_rare.${subtypes[j]}_maf005/CLEAN.rare005.${subtypes[j]}.${pops[i]} \
	      	      --noweb
	done
done




