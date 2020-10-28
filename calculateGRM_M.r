tempDir = "/scratch/genevol/users/lucas/"

listGenes = readRDS ( paste0 ( tempDir , "snpSet_22.rds" ) )

altReadSnps = gaston::read.vcf( paste0 ( tempDir , "chr22.vcf.gz" ) )

#grm_matrix = gaston::as.matrix ( altReadSnps )
#grm_scaled = scale( grm_matrix , center = T , scale = T )

#saveRDS ( grm_scaled , paste0 ( tempDir , "scaledMatrixBk_M22.rds" ) )
grm_scaled = readRDS (paste0 (tempDir , "scaledMatrixBk_M22.rds" ))
matrixDf = as.data.frame ( grm_scaled )
lengthBfRm = ncol(matrixDf)
filteredMatrixDf = matrixDf[ , apply(matrixDf, 2, function(x) !any(is.na(x)))]
matrixAfRm = ncol ( filteredMatrixDf )
rawMatrix = matrix ( t ( filteredMatrixDf [ 1:nrow ( filteredMatrixDf ) , ] ) , nrow = nrow ( filteredMatrixDf ) )

manual_GRM = ( 1 / ncol ( filteredMatrixDf ) ) * ( rawMatrix %*% t ( rawMatrix ) )

saveRDS ( manual_GRM , paste0 ( tempDir , "manualGRM_M22.rds" )  )
saveRDS ( c (lengthBfRm , matrixAfRm) , paste0 ( tempDir , "comparison_M22.rds" ) )
