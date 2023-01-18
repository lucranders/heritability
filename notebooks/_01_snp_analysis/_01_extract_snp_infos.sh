for chr_ in $(seq 1 22); do bcftools query -f '%CHROM\t%POS\t%ID\n' /raid/genevol/vcf_1000G/phase3_20130502_grch38positions/ALL.chr"$chr_"_GRCh38.genotypes.20170504.vcf.gz > test$chr_.txt; done;
