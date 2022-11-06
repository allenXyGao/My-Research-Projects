#!/bin/sh

input=/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output
output=/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/cleaned_dataset
removed_sample_SNPs_dir=/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/removed_Samples_SNPs


module use /projects/rpci/modules
module load plink/1.07

pop="ALLcasecontrol_cohort1 ALLcasecontrol_cohort2 AMLMDScasecontrol_cohort1 AMLMDScasecontrol_cohort2"
for i in $pop
do 
	plink --bfile $input/$i --remove $removed_sample_SNPs_dir/removed.samples_$i \
	      --geno 0.02 --hwe 1e-6 --maf 0.00001 \
              --make-bed \
              --out $output/CLEAN.$i \
	      --noweb	
done


# Q: remove common variants?

# create a new sh.
# read cleaned datset 
# --max-maf 0.01
# out 

# run skato {# deal with missing -> use non-NA}




