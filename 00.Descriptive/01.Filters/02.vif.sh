#!/bin/bash

#PBS -l nodes=1:ppn=4
#PBS -l mem=64gb
#PBS -l walltime=720:00:00
#PBS -q long
#PBS -t 1-22
#PBS -N VifFilters
#PBS -j oe
#PBS -o VifFilters.txt

chr=$PBS_ARRAYID
pathSave=/scratch/genevol/users/lucas/

for vif_ in 1.11 2.5 5; do
for hwe_ in 1e-08 1e-07 1e-06 1e-05 1e-04 0.001 0.01 0.05; do
for maf_ in 0.001 0.005 0.01 0.02 0.03 0.04 0.05; do

input=/raid/genevol/vcf_1000G/phase3_20130502_grch38positions/ALL.chr"$chr"_GRCh38.genotypes.20170504.vcf.gz

outputVif="$pathSave"Vif/Vif_NAfr_chr"$chr"_maf_"$maf_"_hwe_"$hwe_"_vif_"$vif_"
outputVif2="$pathSave"Vif/Vif_Geuvadis_chr"$chr"_maf_"$maf_"_hwe_"$hwe_"_vif_"$vif_"
outputVif3="$pathSave"Vif/Vif_Full_chr"$chr"_maf_"$maf_"_hwe_"$hwe_"_vif_"$vif_"


samp=/raid/genevol/users/lucas/heritability/02.GCTA/data/nAfr.txt
samp2=/raid/genevol/users/lucas/heritability/02.GCTA/data/samples.txt
file_=filter_sample_NAfr_maf_"$maf_"_hwe_"$hwe_"
fileprune="$pathSave"ListSnpFilters/"$file_".txt
file2_=filter_sample_Geuvadis_maf_"$maf_"_hwe_"$hwe_"
fileprune2="$pathSave"ListSnpFilters/"$file2_".txt
file3_=filter_sample_Full_maf_"$maf_"_hwe_"$hwe_"
fileprune3="$pathSave"ListSnpFilters/"$file3_".txt



plink --vcf $input --indep 50 5 $vif_ --extract $fileprune --vcf-half-call missing --keep $samp --out $outputVif &
plink --vcf $input --indep 50 5 $vif_ --extract $fileprune2 --vcf-half-call missing --keep $samp2 --out $outputVif2 &
plink --vcf $input --indep 50 5 $vif_ --extract $fileprune3 --vcf-half-call missing --out $outputVif3

done
done
done