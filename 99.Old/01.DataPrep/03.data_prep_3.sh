# !/bin/bash

# PBS -l nodes=1:ppn=1
# PBS -l mem=16gb
# PBS -l walltime=24:00:00
# PBS -q short
# PBS -N bed2gds
# PBS -t 1-22
# PBS -j oe
# PBS -o log.txt
chr=$PBS_ARRAYID
for chr in $(seq 1 22); do
for version in $(seq 1 3); do

echo ${chr}
echo ${version}

Rscript /raid/genevol/users/lucas/heritability/firstSteps/03.r_bed2gds.r ${chr} ${version}

done
done
