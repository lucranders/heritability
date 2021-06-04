chr_ = commandArgs(TRUE)[1]
version_ =  commandArgs(TRUE)[2]

bed.fn <- paste0("/scratch/genevol/users/lucas/V",version_,"/chr",chr_,".bed")
bim.fn <- paste0("/scratch/genevol/users/lucas/V",version_,"/chr",chr_,".bim")
fam.fn <- paste0("/scratch/genevol/users/lucas/V",version_,"/chr",chr_,".fam")
gds_file <- paste0("/scratch/genevol/users/lucas/V",version_,"/chr",chr_,".gds")

SNPRelate::snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, gds_file)
