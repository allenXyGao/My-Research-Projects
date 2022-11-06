script=/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/scripts/SKATBinary_Robust_maf005     
tmpscript=/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/tmpscript/SKATBinary_Robust
Output=/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/SKATO_test_results_SKATBinary_Robust_BG_maf005
module load R/4.0.4

cd $tmpscript

# need to receive input "pop" from terminal line, e.g. ALLcasecontrol_cohort1
pop=$1
info=/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/SSD_rare_maf005/$pop/File.Info
# number of sets
nSets=$(grep 'TotalNumberOfSets' $info | cut -f1)

start=1

while [ $start -le $nSets ]
do 
	end=$((start-1+100))
	if [ $end -gt $nSets ]
	then
		end=$nSets
	fi
	echo -e "$start\t$end"
 	echo "#!/bin/sh" > SKATBinary_Robust.$pop.$start-$end
	echo "#SBATCH --partition=general-compute" >> SKATBinary_Robust.$pop.$start-$end
	echo "#SBATCH --qos=general-compute" >> SKATBinary_Robust.$pop.$start-$end
	echo "#SBATCH --time=2:00:00" >> SKATBinary_Robust.$pop.$start-$end
	echo "#SBATCH --mem=10000" >> SKATBinary_Robust.$pop.$start-$end
	echo "#SBATCH --nodes=1" >> SKATBinary_Robust.$pop.$start-$end
	echo "#SBATCH --ntasks-per-node=1" >> SKATBinary_Robust.$pop.$start-$end
	echo "#SBATCH --job-name='SKATBinary_Robust$start-$end'" >> SKATBinary_Robust.$pop.$start-$end
	echo "#SBATCH --output=$tmpscript/SKATBinary_Robust.$start-$end.out" >> SKATBinary_Robust.$pop.$start-$end
	echo "##SBATCH --requeue" >> SKATBinary_Robust.$pop.$start-$end
	echo "$script/SKATBinary_Robust.R $pop $start $end > $Output/$pop/logFiles/SKATBinary_Robust.$pop.$start-$end.log" >> SKATBinary_Robust.$pop.$start-$end
	sbatch SKATBinary_Robust.$pop.$start-$end
	start=$((end+1))

done
