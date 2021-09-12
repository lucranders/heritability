# !/bin/bash

# PBS -l nodes=1:ppn=1
# PBS -l mem=16gb
# PBS -l walltime=24:00:00
# PBS -q short
# PBS -N GCTA
# PBS -t 1-22
# PBS -j oe
# PBS -o log.txt

dirFiles=/raid/genevol/users/lucas/heritability/02.GCTA/data
for version in $(seq 1 3); do
Rscript  matrixCorrection.r $version grmNAfr
Rscript  matrixCorrection.r $version grmGeuvadis
Rscript  matrixCorrection.r $version grmFull
done