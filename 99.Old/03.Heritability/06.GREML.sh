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

for version_ in $(seq 1 3); do
grmPath=/scratch/genevol/users/lucas/V$version_
for expr_ in HLA-A HLA-B HLA-C HLA-DPA1 HLA-DPB1 HLA-DQA1 HLA-DQB1 HLA-DRA HLA-DRB1; do
results=/raid/genevol/users/lucas/heritability/03.Heritability/Results/V$version_/$expr_
for num_ in 0 1 2; do
echo ${version_}
echo ${expr_}
echo ${num_}
/raid/genevol/users/lucas/gcta/gcta64 --reml --reml-alg 0 --reml-maxit 1000 --grm $grmPath/grmFull --pheno $dirFiles/gene$expr_.txt --covar $dirFiles/envir$num_.txt --reml-est-fix --out $results/Full$expr_$num_
/raid/genevol/users/lucas/gcta/gcta64 --reml --reml-alg 0 --reml-maxit 1000 --grm $grmPath/grmGeuvadis --pheno $dirFiles/gene$expr_.txt --covar $dirFiles/envir$num_.txt --reml-est-fix --out $results/Geuvadis$expr_$num_
/raid/genevol/users/lucas/gcta/gcta64 --reml --reml-alg 0 --reml-maxit 1000 --grm $grmPath/grmNAfr --pheno $dirFiles/geneNAfr$expr_.txt --covar $dirFiles/envir$num_.txt --reml-est-fix --out $results/NAfr$expr_$num_
done
done
done
