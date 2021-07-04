# !/bin/bash

# PBS -l nodes=1:ppn=1
# PBS -l mem=16gb
# PBS -l walltime=24:00:00
# PBS -q short
# PBS -N bcftools-filter1
# PBS -t 1-22
# PBS -j oe
# PBS -o log.txt


for chr in $(seq 1 22); do

input=/raid/genevol/vcf_1000G/phase3_20130502_grch38positions/ALL.chr"$chr"_GRCh38.genotypes.20170504.vcf.gz

#output=/scratch/genevol/users/lucas/Hardy_NAfr_chr"$chr"
#output2=/scratch/genevol/users/lucas/Hardy_Geuvadis_chr"$chr"
#output3=/scratch/genevol/users/lucas/Hardy_Full_chr"$chr"

outputMaf=/scratch/genevol/users/lucas/Maf/Hardy_NAfr_chr"$chr"
outputMaf2=/scratch/genevol/users/lucas/Maf/Hardy_Geuvadis_chr"$chr"
outputMaf3=/scratch/genevol/users/lucas/Maf/Hardy_Full_chr"$chr"


samp=/raid/genevol/users/lucas/heritability/02.GCTA/data/nAfr.txt
samp2=/raid/genevol/users/lucas/heritability/02.GCTA/data/samples.txt

echo ${chr}
# plink --vcf $input --hardy --vcf-half-call missing --keep $samp --out $output
# plink --vcf $input --hardy --vcf-half-call missing --keep $samp2 --out $output2
# plink --vcf $input --hardy --vcf-half-call missing --out $output3

plink --vcf $input --freq --vcf-half-call missing --keep $samp --out $outputMaf
plink --vcf $input --freq --vcf-half-call missing --keep $samp2 --out $outputMaf2
plink --vcf $input --freq --vcf-half-call missing --out $outputMaf3

done