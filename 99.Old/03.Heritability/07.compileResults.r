library(dplyr)
dfResults = NULL
for(version in c(1,2,3)){
    for(gene_ in c("HLA-A", 'HLA-B', 'HLA-C', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DRA', 'HLA-DRB1')){
        path = paste0('/raid/genevol/users/lucas/heritability/03.Heritability/Results/V',version,'/',gene_)
        elements = list.files(path,full.names=T)
        results = elements[grep('.hsq',elements)]
        for(index in 1:length(results)){
            result = readLines(results[index])
            heritLine = result[5]
            resultP = result[10]
            heritEstimate = as.numeric(strsplit(heritLine,'\t')[[1]][2])
            heritSE = as.numeric(strsplit(heritLine,'\t')[[1]][3])
            heritP = as.numeric(strsplit(resultP,'\t')[[1]][2])
            compResults = as.data.frame(cbind(heritEstimate,heritSE,heritP))
            compResults$geneExp = gene_
            compResults$subPop = ifelse(grepl('Geuvadis',results[index]) , 'Geuvadis' , ifelse(grepl('Full',results[index]) , "Full" , "NAfr"))
            compResults$covs = ifelse(grepl('1\\.hsq',results[index]) , 'Lab+Sex' , 'Lab')
            compResults$VIF = ifelse(grepl('V1',results[index]) , '1.11' , ifelse(grepl('V2',results[index]) , '2.5' , '5'))
            dfResults = rbind(dfResults , compResults)
        }
    }
}
