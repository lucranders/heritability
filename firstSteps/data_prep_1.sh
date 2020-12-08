#!/bin/bash

#PBS -l nodes=1:ppn=1
#PBS -l mem=16gb
#PBS -l walltime=24:00:00
#PBS -q short
#PBS -N bcftools-filter
#PBS -t 1-22
#PBS -j oe
#PBS -o log.txt

chr=$PBS_ARRAYID

input=/raid/genevol/heritability/genotypes_1000g/ALL.chr${chr}_GRCh38.genotypes.20170504.vcf.gz
output=/scratch/genevol/users/lucas/chr${chr}.vcf.gz
samp=/raid/genevol/heritability/samples.txt

bcftools view --samples-file $samp --min-ac=1:minor -o $output $input

Rscript /raid/genevol/users/lucas/heritability/firstSteps/data_prep_1.r ${chr}
