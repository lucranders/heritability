pastaLeitura = "/scratch/genevol/users/lucas/"
pastaSave = "/scratch/genevol/users/lucas/"
arg = commandArgs(trailingOnly=TRUE)

vcf_file <- paste0(pastaLeitura , "chr" , arg ,".vcf.gz")
gds_file <- paste0(pastaSave , "chr" , arg ,".gds")

SeqArray::seqVCF2GDS(vcf_file, gds_file, fmt.import = "GT", verbose = FALSE)

gds <- SeqArray::seqOpen(gds_file)

set.seed(100)
snpset <- SNPRelate::snpgdsLDpruning( gds, method = "corr", slide.max.bp = 10e6, ld.threshold = sqrt(0.1))
saveRDS(snpset, paste0(pastaSave , "snpSet_" , arg , ".rds"))
