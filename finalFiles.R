library(dplyr)

# read and write directory
dirRW = "/scratch/genevol/users/lucas/"

# convert vcf with all chromossomos to gds
vcf_file = paste0( dirRW , "allChr.vcf.gz" )
gds_file = paste0 ( dirRW , "allChr.gds" )
SeqArray::seqVCF2GDS(vcf_file, gds_file, fmt.import = "GT", verbose = FALSE)

# list all files in directory
fileList = list.files ( dirRW , full.names = T ) %>% as.data.frame() %>% rename('fName' = '.')

# filter only files with snps pruned
mainFilesL = fileList %>% filter ( grepl ( "\\.rds" , fName ) )
fullPrunedList = NULL

for ( rdsFile in mainFilesL ){

	fullPrunedList = c ( fullPrunedList , unlist ( rdsFile ) )

}

# save R object containing all snps of interest
saveRDS ( fullPrunedList , paste0 ( dirRW , "fullPrunedList.rds" ) )

# list files to remove
mainFilesR = fileList %>% filter ( !grepl ( "all" , fName ) & ( grepl ( "\\.vcf\\.gz" , fName ) | grepl ( "\\.gds" , fName ) ) )

# build command to remove all separated chromossomes vcf files
rmAll =  "rm"
for (element in mainFilesR$fName){
  
  rmAll = paste ( rmAll , element )
  
}

print ( rmAll )
# effectively remove files
# system ( rmAll )

