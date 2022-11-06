#input=/projects/rpci/qzhu/eschille/BreastGWAS/No_Y_Chrom/maf1.noHM/Input
#output=/projects/rpci/qzhu/eschille/BreastGWAS/No_Y_Chrom/maf1.noHM/Output
input=/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/LDprune_results
output=/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/LDprune_results/LDprune
module load plink/1.90

cd $input

pop="ALLcasecontrol_cohort1 ALLcasecontrol_cohort2 AMLMDScasecontrol_cohort1 AMLMDScasecontrol_cohort2"
for i in $pop
do
	## combine all LDprune result from each chromosome of 1:22, generate the data with all pruned-in SNPs 
	if [ -e $output/LDprune.autosome.maf1.$i.prune.in ]
	then   
        	rm -f $output/LDprune.autosome.maf1.$i.prune.in
	fi
	cat $output/LDprune.autosome.maf1.$i.*.prune.in > $output/LDprune.autosome.maf1.$i.prune.in
	rm $output/LDprune.autosome.maf1.$i.*.prune.*

	plink --bfile $input/maf1.$i.v2 --extract $output/LDprune.autosome.maf1.$i.prune.in --geno 0.02 --maf 0.01 --hwe 1e-6 --make-bed --out $input/$i.maf1.v3.LDprune.autosome.prune.in.$i --noweb
done
