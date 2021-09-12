# !/bin/bash

# PBS -l nodes=1:ppn=1
# PBS -l mem=16gb
# PBS -l walltime=24:00:00
# PBS -q short
# PBS -N bcftools-filter3
# PBS -t 1-22
# PBS -j oe
# PBS -o log.txt


for chr in $(seq 1 22); do

input=/raid/genevol/vcf_1000G/phase3_20130502_grch38positions/ALL.chr"$chr"_GRCh38.genotypes.20170504.vcf.gz
output=/scratch/genevol/users/lucas/chr"$chr"_3.vcf.gz
output2=/scratch/genevol/users/lucas/chr"$chr"_2_3
samp=/raid/genevol/heritability/samples.txt
echo ${chr}
bcftools view --samples-file $samp -o $output $input
plink --vcf $output --hwe 0.0001 --geno 0 --indep 50 5 5 --maf 0.01 --vcf-half-call missing  --out $output2

done
