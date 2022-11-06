i=$1

chip2pca=/projects/rpci/qzhu/qzhu/software/CHIP2PCA/chip2pca.pl
opt=/projects/rpci/qzhu/qzhu/software/CHIP2PCA/current.opt
basename=LDprune_autosome
highLD=/projects/rpci/qzhu/qzhu/software/CHIP2PCA/highLD37.hild
input=/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/EigenStrat_Results/$i/LDprune_autosome-pca
output=/projects/rpci/qzhu/xinyu/mywork/skatO_case_control/output/EigenStrat_Results.iter2/$i

module load plink/1.90
module load R/3.1.2    ## need in chip2pca.pl script.

PATH=$PATH:/projects/rpci/qzhu/qzhu/software/EIGENSOFT/EIG4.2/bin

if [ -e $output ]
then
	rm -fr $output
fi
mkdir -p $output/$basename-cr

cd $output
cp $opt $basename.opt
## only use chr1:22, remove variants in highLD regions, variants with MAF<=1%
perl -ane 'next if(/^#/); $i++; chomp($_); print "$_\tR$i\n"' $highLD > remove.regions

plink --bfile $input/LDprune_autosome-pca2 --exclude remove.regions --range --geno 0.02 --maf 0.01 --hwe 1e-6 --make-bed --out $output/$basename-cr/$basename-cr --noweb

plink --bfile $output/$basename-cr/$basename-cr --recode --out $output/$basename-cr/$basename-cr --noweb

perl $chip2pca $basename snppca

rm -f remove.regions

