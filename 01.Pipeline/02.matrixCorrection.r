path_ =  commandArgs(TRUE)[1]
# Read estimated GRM elements
ReadGRMBin=function(prefix, AllN=F, size=4){
  sum_i=function(i){
    return(sum(1:i))
  }
  BinFileName=paste(prefix,".grm.bin",sep="")
  NFileName=paste(prefix,".grm.N.bin",sep="")
  IDFileName=paste(prefix,".grm.id",sep="")
  id = read.table(IDFileName)
  n=dim(id)[1]
  BinFile=file(BinFileName, "rb");
  grm=readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
  NFile=file(NFileName, "rb");
  if(AllN==T){
    N=readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
  }
  else N=readBin(NFile, n=1, what=numeric(0), size=size)
  i=sapply(1:n, sum_i)
  return(list(diag=grm[i], off=grm[-i], id=id, N=N))
}
# Write elements in binary files
writeGRMBin=function(element_,prefix, size=4){

  BinFileName=paste(prefix,".grm.bin",sep="")
  BinFile=file(BinFileName, "wb");
  grm=writeBin(object = element_,BinFileName, size=size)

}

# Correction for non invertible GRM
fixNonPositiveDefiniteMatrix <- function(origMat) {
    cholStatus <- try(u <- chol(origMat), silent = TRUE)
    cholError <- ifelse(class(cholStatus) == "try-error", TRUE, FALSE)
    newMat <- origMat

    iter <- 0
    while (cholError) {

        iter <- iter + 1
        cat("iteration ", iter, "\n")

        # replace -ve eigen values with small +ve number
        newEig <- eigen(newMat)
        newEig2 <- ifelse(newEig$values < 0, 1e-10, newEig$values)

        # create modified matrix eqn 5 from Brissette et al 2007,
		# inv = transp for eig vectors
        newMat <- newEig$vectors %*% diag(newEig2) %*% t(newEig$vectors)

        # normalize modified matrix eqn 6 from Brissette et al 2007
        newMat <- newMat #/sqrt(diag(newMat) %*% t(diag(newMat)))

        # try chol again
        cholStatus <- try(u <- chol(newMat), silent = TRUE)
        cholError <- ifelse(class(cholStatus) == "try-error", TRUE, FALSE)
    }
    newMat
}

# randomCor, logico, se TRUE, cria as covariancias
# se FALSE, covariancias sao nulas.
blockDiagonalPDMatrix <- function(randomCor,...) {
  matrixList<-list(...)
  if(is.list(matrixList[[1]])) matrixList<-matrixList[[1]]
 
  dimensions<-sapply(matrixList,FUN=function(x) dim(x)[1])
  finalDimension<-sum(dimensions)

  finalMatrix<-matrix(0,nrow=finalDimension,ncol=finalDimension)
  if (randomCor == TRUE)
    finalMatrix <- rcorrmatrix(finalDimension, 1)

  index<-1
  for(k in 1:length(dimensions)){
    finalMatrix[index:(index+dimensions[k]-1),index:(index+dimensions[k]-1)]<-matrixList[[k]]
    index<-index+dimensions[k]
  }
  fixNonPositiveDefiniteMatrix(finalMatrix)
}

options(digits=6)

prefix_ = paste0(path_,'/GCTA')
GRM = ReadGRMBin(prefix_)
M <- matrix(0, dim(GRM[[3]])[1], dim(GRM[[3]])[1])
M[upper.tri(M, diag = FALSE)] = GRM[[2]]
M = t(M)
M[upper.tri(M, diag = FALSE)] = GRM[[2]]
diag(M) = GRM[[1]]
(M[1:4,1:4])

# Matriz corrigida
matrizCorrigida = suppressWarnings(fixNonPositiveDefiniteMatrix(M))
rownames(matrizCorrigida) = GRM[[3]][,1]
colnames(matrizCorrigida) = GRM[[3]][,1]
# matrizCorrigida1 = t(t(matrizCorrigida0)/sqrt(diag(matrizCorrigida0)))
# matrizCorrigida = matrizCorrigida1/sqrt(diag(matrizCorrigida0))
(matrizCorrigida[1:4,1:4])
# vectNew = matrizCorrigida[upper.tri(matrizCorrigida, diag = FALSE)]
# diagNew = diag(matrizCorrigida)
elements = matrizCorrigida[upper.tri(matrizCorrigida, diag = TRUE)]
writeGRMBin(element_ = elements,prefix_)
saveRDS(object = matrizCorrigida,file = paste0(prefix_,'.rds'))
write.table('finished',paste0(path_,'/statusGRMCorrection.txt'))