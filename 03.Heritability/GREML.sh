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

echo ${version}
cd /scratch/genevol/users/lucas/V$version/
/raid/genevol/users/lucas/gcta/gcta64 --mbfile $dirFiles/chrs.txt --keep $dirFiles/nAfr.txt --make-grm --out grmNAfr
/raid/genevol/users/lucas/gcta/gcta64 --mbfile $dirFiles/chrs.txt --keep $dirFiles/samples.txt --make-grm --out grmGeuvadis
/raid/genevol/users/lucas/gcta/gcta64 --mbfile $dirFiles/chrs.txt --make-grm --out grmFull

done
