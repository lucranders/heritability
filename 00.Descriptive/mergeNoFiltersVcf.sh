#!/bin/bash

#PBS -l nodes=1:ppn=1
#PBS -l mem=16gb
#PBS -l walltime=24:00:00
#PBS -q short
#PBS -N MergeNoFilters
#PBS -j oe
#PBS -o MergeNoFilters.txt

input=/raid/genevol/vcf_1000G/phase3_20130502_grch38positions/*.vcf.gz
output=/scratch/genevol/users/lucas/allChrNoFilt.vcf.gz
samp=/raid/genevol/heritability/samples.txt

bcftools concat -o $output $input