#!/bin/bash

#PBS -l nodes=1:ppn=1
#PBS -l mem=16gb
#PBS -l walltime=24:00:00
#PBS -q short
#PBS -N bcftools-filter
#PBS -t 1-19
#PBS -j oe
#PBS -o log.txt

chr=$PBS_ARRAYID

Rscript /raid/genevol/users/lucas/heritability/manualDistMPartsSerial.r ${chr}
