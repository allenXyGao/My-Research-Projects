## generate a text file containing samples, reported race, eigen values of each PC axis from the output of eigenstrat:

program=/projects/rpci/qzhu/qliu7/myscript/CHIP2PCA.analysis.R
dir=/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/EigenStrat_Results
basename=LDprune_autosome
module load R/3.1.2

pop="ALLcasecontrol_cohort1 ALLcasecontrol_cohort2 AMLMDScasecontrol_cohort1 AMLMDScasecontrol_cohort2"
for i in $pop
do
	path=$dir/$i/$basename-pca
	for j in {1..2}
	do 
		cut -d" " -f1,2 $path/$basename-pca$j.fam > $path/pca$j
		awk -v p=$i '{print $1 ":" $2 "\t" p}' $path/pca$j > $path/pca$j.eth
		R --slave --no-save --args input="$path/$basename-pca$j.evec" population="$path/pca$j.eth" output="$dir/$i/$basename-pca$j.eth.out" < $program
	done
	# rm $path/pca*
done

#1. write a loop for each population.
#2. generate eth/race info for each population.  
#3. use pca1 and pca2 separately, generate 2 output files. 
# 4. mark removed samples from pca2, and mark them in the scatterplot of pca1.(in scatterplot step)

