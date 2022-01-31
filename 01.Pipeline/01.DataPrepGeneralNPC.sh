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


input=/raid/genevol/vcf_1000G/phase3_20130502_grch38positions/ALL.chr"$chr"_GRCh38.genotypes.20170504.vcf.gz

bed=$tempPath/chr"$chr"
$query_