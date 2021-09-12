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

descriptiveAnalisysGrm = function(version_,matrixName_,list){
    
    prefix = paste0('/scratch/genevol/users/lucas/V' , version_ , '/' , matrixName_)

    itensAfr = ReadGRMBin(prefix)
    M <- matrix(0, dim(itensAfr[[3]])[1], dim(itensAfr[[3]])[1])
    M[upper.tri(M, diag = FALSE)] = itensAfr[[2]]
    M = t(M)
    M[upper.tri(M, diag = FALSE)] = itensAfr[[2]]
    diag(M) = itensAfr[[1]]
    rownames(M) = itensAfr[[3]][,1]
    colnames(M) = itensAfr[[3]][,1]
    (M[1:4,1:4])

    GRM = M[list,list]
    return(GRM)

}


eigenDecomp = function(GRM,label_,pops_){
    
    eig = eigen(GRM)
    
    eigenValuesDf = cbind(eig$values,1:length(eig$values) ) %>% as.data.frame()
    colnames ( eigenValuesDf ) = c ( "eigen" , "order" )
    eigenValuesDf$cummulativeEigen = cumsum( eigenValuesDf$eigen ) / sum(eigenValuesDf$eigen)
    eigenValuesDf$flag = ifelse ( eigenValuesDf$eigen > mean(eigenValuesDf$eigen) , '>1' , '<=1' )
    eigenValuesDf$flag2 = ifelse ( eigenValuesDf$cummulativeEigen > .5 , '>1/2' , '<=1/2' )
    eigenValuesDf$label = label_


    pcadf <- as.data.frame(eig$vectors[,1:2]) %>%
    as_tibble() %>%
    mutate(sampleid = rownames(GRM) , label = label_) %>%
    inner_join(pops_, ., by = "sampleid")


    return(list(eigenValuesDf,pcadf))

}

listGeuvadis = read.table('/raid/genevol/heritability/samples.txt',header = F)
colnames(listGeuvadis) = "samples"
# informacao de expressao
hlaExp = readr::read_tsv("/raid/genevol/heritability/hla_expression.tsv") %>% rename(sampleid = subject_id) %>% filter(sampleid %in% listGeuvadis$samples)
pops = hlaExp %>%  distinct(sampleid , pop)
popsNAfr = pops %>% filter(pop != "YRI")


dfEigenVal = NULL
dfEigenVect = NULL
for(version in 1:3){

    assign(paste0("kinshipV",version,"Geuvadis") ,  descriptiveAnalisysGrm(version,'grmGeuvadis',listGeuvadis$samples))
    assign(paste0("kinshipV",version,"Full") ,  descriptiveAnalisysGrm(version,'grmFull',listGeuvadis$samples))
    assign(paste0("kinshipV",version,"NAfr") ,  descriptiveAnalisysGrm(version,'grmNAfr',popsNAfr$sampleid))

    valuesV1Geuvadis = eigenDecomp( get(paste0("kinshipV",version,"Geuvadis")) , paste0('Geuvadis',version) , pops ) 
    valuesV1Full = eigenDecomp( get(paste0("kinshipV",version,"Full") ) , paste0('Full',version) , pops )
    valuesV1NAfr = eigenDecomp( get(paste0("kinshipV",version,"NAfr") ) , paste0('NAfr',version) , pops )

    dfEigenVal = rbind(dfEigenVal , valuesV1Geuvadis[[1]] )
    dfEigenVal = rbind(dfEigenVal , valuesV1Full[[1]] )
    dfEigenVal = rbind(dfEigenVal , valuesV1NAfr[[1]] )
    dfEigenVect = rbind(dfEigenVect , valuesV1Geuvadis[[2]] )
    dfEigenVect = rbind(dfEigenVect , valuesV1Full[[2]] )
    dfEigenVect = rbind(dfEigenVect , valuesV1NAfr[[2]] )

}

p1 = dfEigenVect %>% ggplot(aes(x = V1,y = V2 , colour = pop)) +
geom_point() +
facet_wrap(~label,ncol = 3) +
labs(x = "PC1" , y = "PC2" , colour = "Pop" , title = "Comparison between methods") +
theme_classic()
ggsave('/raid/genevol/users/lucas/heritability/03.Heritability/PCcomparison.png' , p1 , device = 'png')

