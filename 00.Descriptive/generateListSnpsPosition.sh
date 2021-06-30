# !/bin/bash

# PBS -l nodes=1:ppn=1
# PBS -l mem=16gb
# PBS -l walltime=24:00:00
# PBS -q short
# PBS -N bcftools-references snps
# PBS -t 1-22
# PBS -j oe
# PBS -o log.txt

pathRead=/raid/genevol/vcf_1000G/phase3_20130502_grch38positions
pathSave=/scratch/genevol/users/lucas

for num in $(seq 1 22); do

echo ${num}
bcftools query -f '%CHROM %ID %POS\n' $pathRead/ALL.chr"$num"_GRCh38.genotypes.20170504.vcf.gz -o $pathSave/basicInfoChr$num.txt

done
