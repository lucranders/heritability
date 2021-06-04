tempDir = "/scratch/genevol/users/lucas/"

#listGenes = readRDS ( paste0 ( tempDir , "fullPrunedList.rds" ) )

#altReadSnps = gaston::read.vcf( paste0 ( tempDir , "allChr.vcf.gz" ) )

#grm_matrix = gaston::as.matrix ( altReadSnps )
#grm_scaled = scale( grm_matrix , center = T , scale = T )

grm_scaled = readRDS (paste0 (tempDir , "scaledMatrixBk.rds" ))
matrixDf = as.data.frame ( grm_scaled )
lengthBfRm = ncol(matrixDf)
print("read Df: ok")
filteredMatrixDf = matrixDf[ , apply(matrixDf, 2, function(x) !any(is.na(x)))]
print("filter df: ok")
matrixAfRm = ncol ( filteredMatrixDf )
rawMatrix = matrix ( t ( filteredMatrixDf [ 1:nrow ( filteredMatrixDf ) , ] ) , nrow = nrow ( filteredMatrixDf ) )
print("Df as matrix: ok")
manual_GRM = ( 1 / ncol ( filteredMatrixDf ) ) * ( rawMatrix %*% t ( rawMatrix ) )
print("GRM: ok")

saveRDS ( manual_GRM , paste0 ( tempDir , "manualGRM.rds" )  )
saveRDS ( c (lengthBfRm , matrixAfRm) , paste0 ( tempDir , "comparison.rds" ) )
