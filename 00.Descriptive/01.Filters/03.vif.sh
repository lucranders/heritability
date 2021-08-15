#!/bin/bash

#PBS -l nodes=1:ppn=1
#PBS -l mem=16gb
#PBS -l walltime=24:00:00
#PBS -q short
#PBS -N VifFilters
#PBS -j oe
#PBS -o VifFilters.txt

pathRead=/raid/genevol/vcf_1000G/phase3_20130502_grch38positions
pathSave=/scratch/genevol/users/lucas/Vif/

for chr in $(seq 1 22); do
for vif_ in 1.11 2.5 5; do
echo ${chr}
echo ${vif_}

input=/raid/genevol/vcf_1000G/phase3_20130502_grch38positions/ALL.chr"$chr"_GRCh38.genotypes.20170504.vcf.gz

outputVif="$pathSave"Vif_NAfr_chr"$chr"_vif_"$vif_"
outputVif2="$pathSave"Vif_Geuvadis_chr"$chr"_vif_"$vif_"
outputVif3="$pathSave"Vif_Full_chr"$chr"_vif_"$vif_"


samp=/raid/genevol/users/lucas/heritability/02.GCTA/data/nAfr.txt
samp2=/raid/genevol/users/lucas/heritability/02.GCTA/data/samples.txt


plink --vcf $input --indep 50 $vif_ --vcf-half-call missing --keep $samp --out $outputVif
plink --vcf $input --indep 50 $vif_ --vcf-half-call missing --keep $samp2 --out $outputVif2
plink --vcf $input --indep 50 $vif_ --vcf-half-call missing --out $outputVif3

done
done