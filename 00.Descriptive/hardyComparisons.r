suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(ggplot2)))

root = '/scratch/genevol/users/lucas/'

hweNafrDat = NULL
hweGeuvadisDat = NULL
hweFullDat = NULL

for(i in 1:22){
    i=22
    auxNafr = data.table::fread(paste0(root,'Hardy_NAfr_chr',i,'.hwe'), header = T,select = c('CHR' , 'SNP' , 'P'))
    auxGeuvadis = data.table::fread(paste0(root,'Hardy_Geuvadis_chr',i,'.hwe'), header = T,select = c('CHR' , 'SNP' , 'P'))
    auxFull = data.table::fread(paste0(root,'Hardy_Full_chr',i,'.hwe'), header = T,select = c('CHR' , 'SNP' , 'P'))
    auxRef = data.table::fread(paste0(root,'basicInfoChr',i,'.txt'), header = F)
    colnames(auxRef) = c('CHR','SNP','BP') 
    
    hweNafrDat = rbind(hweNafrDat , merge(auxNafr , auxRef))
    hweGeuvadisDat = rbind(hweGeuvadisDat , merge(auxGeuvadis , auxRef))
    hweFullDat = rbind(hweFullDat , merge(auxFull , auxRef))

}
