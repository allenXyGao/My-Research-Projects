#input=/projects/rpci/qzhu/eschille/BreastGWAS/No_Y_Chrom/maf1.noHM/Input
#output=/projects/rpci/qzhu/eschille/BreastGWAS/No_Y_Chrom/maf1.noHM/Output
input=/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output
output=/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/LDprune_results
module load plink/1.90

cd $output

pop="ALLcasecontrol_cohort1 ALLcasecontrol_cohort2 AMLMDScasecontrol_cohort1 AMLMDScasecontrol_cohort2"
for i in $pop
do
	#plink --bfile $input/$i --remove $output/Sexcheck/sexcheck.x.problem.$i --geno 0.02 --maf 0.01 --hwe 1e-6 --make-bed --out $input/maf1.$i.v2 --noweb
	plink --bfile $input/$i --geno 0.02 --maf 0.01 --hwe 1e-6 --make-bed --out $output/maf1.$i.v2 --noweb	

	for j in {1..22}
	do

		#plink --bfile $output/maf1.$i.v2 --chr $j --indep-pairwise 1500 150 0.3 --out /projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/LDprune_results/LDprune/LDprune.autosome.maf1.$i.$j --noweb
		plink --bfile $output/maf1.$i.v2 --chr $j --indep-pairwise 1500 150 0.3 --out $output/LDprune/LDprune.autosome.maf1.$i.$j --noweb

	done
done
