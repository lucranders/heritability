# !/bin/bash

# PBS -l nodes=1:ppn=1
# PBS -l mem=16gb
# PBS -l walltime=24:00:00
# PBS -q short
# PBS -N vcf2bed
# PBS -t 1-22
# PBS -j oe
# PBS -o log.txt

for chr in $(seq 1 22); do
for version in $(seq 1 3); do

input=/raid/genevol/vcf_1000G/phase3_20130502_grch38positions/ALL.chr"$chr"_GRCh38.genotypes.20170504.vcf.gz
fileprune=/scratch/genevol/users/lucas/chr"$chr"_2_"$version".prune.in
pruneddata=/scratch/genevol/users/lucas/V"$version"/chr"$chr"
echo ${chr}
echo ${fileprune}
echo ${version}
echo ${pruneddata}
plink --vcf $input --vcf-half-call missing --extract $fileprune --make-bed --out $pruneddata

done
done
