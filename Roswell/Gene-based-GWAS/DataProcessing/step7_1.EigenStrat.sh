chip2pca=/projects/rpci/qzhu/qzhu/software/CHIP2PCA/chip2pca.pl
opt=/projects/rpci/qzhu/qzhu/software/CHIP2PCA/current.opt
input=/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/LDprune_results
basename=LDprune_autosome
highLD=/projects/rpci/qzhu/qzhu/software/CHIP2PCA/highLD37.hild
i=$1
output=/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/EigenStrat_Results/$i


#$i.maf1.v3.LDprune.autosome.prune.in.$i


module load plink/1.90
module load R/3.3.2    ## need in chip2pca.pl script.export 

PATH=$PATH:/projects/rpci/qzhu/qzhu/software/EIGENSOFT/EIG-6.1.4/bin

if [ -e $output ]
then
	rm -fr $output
fi
mkdir -p $output/$basename-cr

cd $output
cp $opt $basename.opt 
## only use chr1:22, remove variants in highLD regions, variants with MAF<=1%
perl -ane 'next if(/^#/); $i++; chomp($_); print "$_\tR$i\n"' $highLD > remove.regions

## prepare binary files in format that chip2pca required. 
## 1) copy .bed file
cp $input/$i.maf1.v3.LDprune.autosome.prune.in.$i.bed $input/$i.maf1.v3.LDprune.autosome.prune.in.$i.chip2pca.bed

## 2) correct .bim file: the col3 (genetic distance in Morgan) to 0. force the chip2pca to not use this value to do the calculation. 
perl -ane 'print "$F[0]\t$F[1]\t0\t$F[3]\t$F[4]\t$F[5]\n"' $input/$i.maf1.v3.LDprune.autosome.prune.in.$i.bim > $input/$i.maf1.v3.LDprune.autosome.prune.in.$i.chip2pca.bim

## 3) correct .fam file: cut the col2 of the IID names shorter; correct the col6 from -9 to 1, both means missing phenotype. 
perl -ane 'print "$F[0] $F[1] $F[2] $F[3] $F[4] 1\n"' $input/$i.maf1.v3.LDprune.autosome.prune.in.$i.fam > $input/$i.maf1.v3.LDprune.autosome.prune.in.$i.chip2pca.fam

##pop-specific
# here only exclude remove.regions, no remove.PLINK.v3.$i.txt
plink --bfile $input/$i.maf1.v3.LDprune.autosome.prune.in.$i.chip2pca --exclude remove.regions --range --geno 0.02 --maf 0.01 --hwe 1e-6 --make-bed --out $output/$basename-cr/long-cr --noweb


## shorten IID names
#awk '{gsub(/-/," ",$2)}7' $output/$basename-cr/long-cr.fam > $output/$basename-cr/temp.fam
awk '{gsub(/[A-Z]-/," ",$2)}7' $output/$basename-cr/long-cr.fam > $output/$basename-cr/temp.fam
perl -ane 'print "$F[0] $F[2] $F[3] $F[4] $F[5] $F[6] \n"' $output/$basename-cr/temp.fam > $output/$basename-cr/long-cr.fam

plink --bfile $output/$basename-cr/long-cr --make-bed --out $output/$basename-cr/$basename-cr --noweb
plink --bfile $output/$basename-cr/long-cr --recode --out $output/$basename-cr/$basename-cr --noweb

perl $chip2pca $basename snppca

rm -f remove.regions inbreeding.IBD.remove.$i



# made a mistake here
#awk '{gsub(/[A-Z]-/," ",$2)}7' /projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/LDprune_results/AMLMDScasecontrol_cohort1.maf1.v3.LDprune.autosome.prune.in.AMLMDScasecontrol_cohort1.chip2pca.fam | head










