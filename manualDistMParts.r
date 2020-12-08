library(dplyr)

tempDir = "/scratch/genevol/users/lucas/"
saveDir = "/raid/genevol/users/lucas/heritability/plots/"
arg = commandArgs(trailingOnly=TRUE)

grmScaled = readRDS (paste0(tempDir , "scaledMatrixBk.rds"))

listSnps = readRDS ( paste0 ( tempDir , "snpSet_" , arg , ".rds" ) )

assign ( paste0 ( "matrixPtChr" , arg ) , grmScaled [ , listSnps [[1]]] )

lengthB = length(listSnps[[1]])
noNa = function(x) !any(is.na(x))
naFreeDf = get ( paste0 ( "matrixPtChr" , arg ) ) %>% as.data.frame ( ) %>% select_if ( noNa ) 
naFreeColumns = colnames(naFreeDf)
lengthA = length ( naFreeColumns )

naFreeMatrix = get ( paste0 ( "matrixPtChr" , arg ) ) [ , naFreeColumns ]

DistMatrix = naFreeMatrix %*% t ( naFreeMatrix )
listLengths = list ( lengthB , lengthA )

saveRDS ( DistMatrix , paste0 ( tempDir , "distMatrixPt", arg , ".rds" ) )
saveRDS ( listLengths , paste0 ( tempDir , "listLengths", arg , ".rds" ) )
