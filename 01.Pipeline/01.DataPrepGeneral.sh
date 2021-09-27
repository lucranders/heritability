
#!/bin/bash

#PBS -l nodes=1:ppn=16
#PBS -l mem=64gb
#PBS -l walltime=24:00:00
#PBS -q short
#PBS -t 1-22
#PBS -N bedMaker
#PBS -j oe
#PBS -o bedMaker.txt


chr=$PBS_ARRAYID


filteredSnps=/scratch/genevol/users/lucas/Vif/Vif_"$pop"_chr"$chr"_maf_"$maf_"_hwe_"$hwe_"_vif_"$vif_".prune.in
input=/raid/genevol/vcf_1000G/phase3_20130502_grch38positions/ALL.chr"$chr"_GRCh38.genotypes.20170504.vcf.gz

tempPath=/scratch/genevol/users/lucas/TempBed_pop_"$pop"_maf_"$maf_"_hwe_"$hwe_"_vif_"$vif_"
bed=$tempPath/chr"$chr"

if [ $pop = NAfr ]; then
    filesamp=/raid/genevol/users/lucas/heritability/01.Pipeline/Samples/nAfr.filt
    plink --vcf $input --vcf-half-call missing --keep $filesamp --extract $filteredSnps --make-bed --out $bed
elif [ $pop = Geuvadis ]; then
    filesamp=/raid/genevol/users/lucas/heritability/01.Pipeline/Samples/samples.filt
    plink --vcf $input --vcf-half-call missing --keep $filesamp --extract $filteredSnps --make-bed --out $bed
else
    plink --vcf $input --vcf-half-call missing --extract $filteredSnps --make-bed --out $bed
fi
