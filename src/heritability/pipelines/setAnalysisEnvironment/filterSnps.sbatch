#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --partition=short
#SBATCH --array=1-22
#SBATCH --time=480
#SBATCH --mem=16gb
#SBATCH --job-name=bedMaker
#SBATCH --nodelist=darwin
#SBATCH --output=/scratch/genevol/users/lucas/out/bedMaker.txt


chr=$SLURM_ARRAY_TASK_ID


input=$pathVcf/ALL.chr"$chr"_GRCh38.genotypes.20170504.vcf.gz

bed=$tempPath/chr"$chr"
eval $query_