script=/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/scripts

pop="ALLcasecontrol_cohort1 ALLcasecontrol_cohort2 AMLMDScasecontrol_cohort1 AMLMDScasecontrol_cohort2"
#pop="ALLcasecontrol_cohort1"
for i in $pop
do
        echo "#!/bin/sh" > EigenStrat_on_slurm.$i
        echo "#SBATCH --partition=general-compute" >> EigenStrat_on_slurm.$i
	echo "#SBATCH --mem=128000" >> EigenStrat_on_slurm.$i
        echo "#SBATCH --time=10:00:00" >> EigenStrat_on_slurm.$i
        echo "#SBATCH --nodes=1" >> EigenStrat_on_slurm.$i
        echo "#SBATCH --ntasks-per-node=8" >> EigenStrat_on_slurm.$i
        echo "#SBATCH --job-name='EigenStrat.$i'" >> EigenStrat_on_slurm.$i
        echo "#SBATCH --output=EigenStrat.$i.ccr.out" >> EigenStrat_on_slurm.$i
        echo "#SBATCH --mail-user=xinyugao@buffalo.edu" >> EigenStrat_on_slurm.$i
        echo "#SBATCH --mail-type=ALL" >> EigenStrat_on_slurm.$i
        echo "##SBATCH --requeue" >> EigenStrat_on_slurm.$i
        echo "$script/step7_1.EigenStrat.sh $i" >> EigenStrat_on_slurm.$i
        sbatch EigenStrat_on_slurm.$i
done



