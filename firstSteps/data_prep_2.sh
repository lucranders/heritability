#!/bin/bash

#PBS -l nodes=1:ppn=1
#PBS -l mem=16gb
#PBS -l walltime=24:00:00
#PBS -q short
#PBS -N fileAdjust
#PBS -j oe
#PBS -o fAdjust.txt

input=/scratch/genevol/users/lucas/*.vcf.gz
output=/scratch/genevol/users/lucas/allChr.vcf.gz

bcftools concat -o $output $input

Rscript /raid/genevol/users/lucas/heritability/firstSteps/data_prep_2.r
