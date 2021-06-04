# !/bin/bash

# PBS -l mem=16gb
# PBS -l walltime=24:00:00
# PBS -q short
# PBS -N DataPrep
# PBS -t 1-22
# PBS -j oe
# PBS -o log.txt
dataPrepPath=/raid/genevol/users/lucas/heritability/01.DataPrep
makeGrm=/raid/genevol/users/lucas/heritability/02.GCTA
echo "Part 1 - 1"
sh $dataPrepPath/01.data_prep_1_1.sh &
sh $dataPrepPath/01.data_prep_1_2.sh &
sh $dataPrepPath/01.data_prep_1_3.sh
echo "Part 2"
sh $dataPrepPath/02.data_prep_2.sh
echo "Part 3"
sh $makeGrm/makeGRM.sh