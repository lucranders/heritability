suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(ggplot2)))

rootSave = '/raid/genevol/users/lucas/heritability/00.Descriptive/Images/Phenotype/'
rootSaveTables = '/raid/genevol/users/lucas/heritability/00.Descriptive/Tables/'
listGeuvadis = read.table('/raid/genevol/heritability/samples.txt',header = F)
colnames(listGeuvadis) = "samples"
hlaExp = readr::read_tsv("/raid/genevol/heritability/hla_expression.tsv") %>% rename(sampleid = subject_id) %>% filter(sampleid %in% listGeuvadis$samples)

# List of all unique genes of interest
genesLoop = sort(unique(hlaExp$gene_name))
# Creating data frame to plot assymetry
dfSymmetryPhen = as.data.frame(NULL)
dfGreatesValues = as.data.frame(1:5)
idsHigher = c()

for( idx_ in 1:length(genesLoop) ){
    # For each gene:
    gene_ = genesLoop[idx_]
    # Filter gene associated information
    aux = hlaExp %>% filter(gene_name == gene_)
    # Selecting samples in top 10 tpm
    top10 = aux %>% arrange(-tpm) %>% mutate(n = 1:n()) %>% filter(n <= 10) %>% select(c('sampleid','tpm'))
    idsHigher = c(idsHigher,top10$sampleid)
    colnames(top10) = c(paste0(gene_,'_sample') , paste0(gene_,'_tpm'))
    dfGreatesValues = cbind(dfGreatesValues,top10)
    # Sort in ascending order gene expression (tpm)
    tpmSort = sort(aux$tpm)
    # Separate ordered set in two parts:
    # the 50% lower values - valuesInf
    #  and the 50% greater values - valuesSup
    # Then measure the distance to the median
    lengthSamp = length(tpmSort)
    sizeParts = round(lengthSamp/2)
    valuesInf = quantile(aux$tpm,.5) - tpmSort[1:sizeParts]
    valuesSup = tpmSort[seq(from = lengthSamp , to = (lengthSamp - sizeParts + 1))] - quantile(aux$tpm,.5)
    # create auxiliary df to append info
    dfAppend = as.data.frame(valuesInf)
    dfAppend$Sup_ = valuesSup
    dfAppend$gene_ = gene_
    # Append in final df
    dfSymmetryPhen = rbind(dfSymmetryPhen , dfAppend)
}

# Rename column
colnames(dfSymmetryPhen)[1] = 'Inf_'

# Univariate analysus:
# Saving tables containing top 10 samples by gene
idsHigher = as.data.frame(idsHigher)
idsHigher %>% group_by(idsHigher) %>% summarise(n=n()) %>% ungroup() %>% arrange(-n)
write.table(dfGreatesValues[,2:7],'/raid/genevol/users/lucas/heritability/00.Descriptive/Tables/top10TpmPt1.txt',sep = '&',quote = F , row.names = F,col.names = T)
write.table(dfGreatesValues[,8:13],'/raid/genevol/users/lucas/heritability/00.Descriptive/Tables/top10TpmPt2.txt',sep = '&',quote = F , row.names = F,col.names = T)
write.table(dfGreatesValues[,14:19],'/raid/genevol/users/lucas/heritability/00.Descriptive/Tables/top10TpmPt3.txt',sep = '&',quote = F , row.names = F,col.names = T)

# Plot assymetry
plotSymmetry = dfSymmetryPhen %>% ggplot(aes(x = Inf_ , y = Sup_ , colour = gene_)) +
geom_point(alpha = .8) +
geom_abline(intercept = 0 , slope = 1) +
theme_bw() +
theme( 
      legend.position='none',
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
facet_wrap(~gene_ , scales = 'free')+
labs(x = "Distance from median (lower)" , y= 'Distance from median (upper)' , colour = 'Gene')

# Plot phenotypes by gene (boxplot)
boxplotGenes = hlaExp %>% ggplot(aes(x = gene_name , y=tpm , fill = gene_name)) +
geom_boxplot(alpha = .8) +
theme_bw() +
theme( 
      legend.position='none',
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
labs(x = 'Gene',y= "Transcriptions per million" , fill = 'Gene')

# Bind plots 
bindPlots = gridExtra::grid.arrange(boxplotGenes,plotSymmetry, ncol = 1)
# Save plots
ggsave(paste0(rootSave,'PhenoDensitySym',".png"), bindPlots, bg = "transparent",width=10, height=8, dpi=300)

# Multivariate analysis
# Transforming original dataset in horizontal position to analyse
horizontalExpDf = tidyr::spread(hlaExp,gene_name,tpm)
# Filter columns of interest - gene expressions and samples
horizontalTpm = horizontalExpDf[,grepl('HLA|sample',colnames(horizontalExpDf))]

# Calculate means and covariance matrix
meanExp = colMeans(horizontalTpm[,2:10])
covExp = cov(horizontalTpm[,2:10])
# Calculate Mahalanobis distance
mahalanobisD = mahalanobis(horizontalTpm[,2:10],center = meanExp,cov = covExp)
horizontalExpDf$Mahalanobis = mahalanobisD
horizontalExpDf$idx = 1:445
# Defining threshold on 99% chi-squared (9 df) quantile 
cutMahalanobis = qchisq(.99, df = 9)
# Create dummy variable to define colours in the plot
horizontalExpDf = horizontalExpDf %>% mutate(threshold = ifelse(Mahalanobis > cutMahalanobis , 'yes' , 'no'))

dfOrderedMahalanobis = as.data.frame( sort ( horizontalExpDf$Mahalanobis) )
dfOrderedMahalanobis$idx = ((1:nrow(dfOrderedMahalanobis)) - .5)/nrow(dfOrderedMahalanobis)
dfOrderedMahalanobis$chisqQ = qchisq(dfOrderedMahalanobis$idx, df = 9)
dfOrderedMahalanobisPlt = dfOrderedMahalanobis[,c(1,3)]
colnames(dfOrderedMahalanobisPlt) = c('Observed','Expected')
dfOrderedMahalanobisPlt = dfOrderedMahalanobisPlt %>% mutate(Threshold = ifelse(Observed > cutMahalanobis , 'yes' , 'no'))

confidenceBandsMahalanobis = function(n,qtTrials = 1000,idxInf = 25 , idxSup = 975){

    sample_ = matrix(0,n,qtTrials)
    limInf = numeric(n)
    limSup = numeric(n)
    #
    set.seed(9297791)
    for(i in 1:qtTrials){
        sample_[,i] = rchisq(n,df = 9)
        sample_[,i] = sort(sample_[,i]) 
        }
    #
    for(i in 1:n){
        quantiles_ = sort(sample_[i,])
        limInf[i]  = quantiles_[idxInf]
        limSup[i]  = quantiles_[idxSup]
        
        }
    return(list(limInf,limSup))
}

intervalFull = confidenceBandsMahalanobis(n = nrow(dfOrderedMahalanobisPlt) , idxInf = 50 , idxSup = 9950 , qtTrials = 10000)
dfOrderedMahalanobisPlt$LimInf = intervalFull[[1]]
dfOrderedMahalanobisPlt$LimSup = intervalFull[[2]]


plotQQFull = dfOrderedMahalanobisPlt %>% 
ggplot() +
geom_point( aes ( x = Expected , y = Observed , colour = Threshold) ) +
geom_line ( aes ( x = Expected , y = LimInf ) , linetype = 2 , colour = 'red' ) +
geom_line ( aes ( x = Expected , y = LimSup ) , linetype = 2 , colour = 'red' ) +
geom_abline( slope = 1 , intercept = 0 ) +
theme_bw() +
theme( 
      legend.position='none',
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) 

ggsave(paste0(rootSave,'mahalanobisQq',".png"), plotQQFull, bg = "transparent",width=10, height=8, dpi=300)


filterOut = horizontalExpDf %>% filter(threshold == 'no')
newMaha = filterOut[grepl('HLA|sample',colnames(horizontalExpDf))]

meanExpF = colMeans(newMaha[,2:10])
covExpF = cov(newMaha[,2:10])
# Calculate Mahalanobis distance
mahalanobisDF = mahalanobis(newMaha[,2:10],center = meanExpF,cov = covExpF)
filterOut$Mahalanobis = mahalanobisDF

dfOrderedMahalanobisF = as.data.frame( sort ( filterOut$Mahalanobis) )
dfOrderedMahalanobisF$idx = ((1:nrow(dfOrderedMahalanobisF)) - .5)/nrow(dfOrderedMahalanobisF)
dfOrderedMahalanobisF$chisqQ = qchisq(dfOrderedMahalanobisF$idx, df = 9)
dfOrderedMahalanobisPltF = dfOrderedMahalanobisF[,c(1,3)]
colnames(dfOrderedMahalanobisPltF) = c('Observed','Expected')

intervalFilt = confidenceBandsMahalanobis(n = nrow(dfOrderedMahalanobisPltF)  , idxInf = 50 , idxSup = 9950 , qtTrials = 10000)
dfOrderedMahalanobisPltF$LimInf = intervalFilt[[1]]
dfOrderedMahalanobisPltF$LimSup = intervalFilt[[2]]


plotQQFilt =  dfOrderedMahalanobisPltF %>% 
ggplot() +
geom_point( aes ( x = Expected , y = Observed ) ) +
geom_line ( aes ( x = Expected , y = LimInf ) , linetype = 2 , colour = 'red' ) +
geom_line ( aes ( x = Expected , y = LimSup ) , linetype = 2 , colour = 'red' ) +
geom_abline( slope = 1 , intercept = 0 ) +
theme_bw() +
theme( 
      legend.position='none',
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) 

ggsave(paste0(rootSave,'mahalanobisQqFilt',".png"), plotQQFilt, bg = "transparent",width=10, height=8, dpi=300)


# Plot Mahalanobis Distance by index
plotMahalanobis = horizontalExpDf %>% 
ggplot(aes(x = idx , y = Mahalanobis , colour =  threshold)) +
geom_point() +
geom_hline(yintercept=cutMahalanobis) +
ggrepel::geom_text_repel(data=subset(horizontalExpDf, Mahalanobis > cutMahalanobis),
            aes(idx,Mahalanobis,label=sampleid)) +
theme_bw() +
theme( 
      legend.position='none',
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )

ggsave(paste0(rootSave,'mahalanobis',".png"), plotMahalanobis, bg = "transparent",width=10, height=8, dpi=300)

# Standardize gene expressions
horizontalTpm[,2:10] = scale(horizontalTpm[,2:10],center = T , scale = T)
# Create biplot - PCA for XX'and X'X
biplotFit = princomp(x = horizontalTpm[,2:10], cor=TRUE)
biplotDf = data.frame(sampleid = horizontalTpm[,1] , biplotFit$scores)
# Merge df with dummy variable based on Mahalanobis distance and threshold
biplotDf = merge(biplotDf,horizontalExpDf[,c('sampleid','threshold')])
pcLoads = data.frame(biplotFit$loadings[,1:9])
scaleBip = 2

# Plot biplot
plotBiplot = biplotDf %>% ggplot() + 
geom_text(data=subset(biplotDf, threshold == 'yes'),
            aes(Comp.1,Comp.2,label=sampleid,colour = threshold),alpha = .7) +
geom_text(data=subset(biplotDf, threshold == 'no'),
            aes(Comp.1,Comp.2,label=sampleid,colour = threshold),alpha = .5) +
geom_hline(yintercept = 0, size=.2) + 
geom_vline(xintercept = 0, size=.2) +
theme_bw() +
theme( 
      legend.position='none',
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
geom_segment(data = pcLoads , aes(x=0, y=0,xend = scaleBip*Comp.1 , yend = scaleBip*Comp.2), arrow=arrow(length=unit(0.1,"cm")), alpha=0.75) +
ggrepel::geom_text_repel(data = pcLoads , aes(x=scaleBip*Comp.1, y=scaleBip*Comp.2, label=row.names(pcLoads) , size = 5, vjust=1) , ) +
labs(x = 'Principal Component 1' , y = 'Principal Component 2')


ggsave(paste0(rootSave,'biplotGeneExp',".png"), plotBiplot, bg = "transparent",width=10, height=8, dpi=300)

# Selecting 'outliers'
samplesDf = horizontalExpDf %>% filter(threshold == 'yes') %>% select(sampleid,Mahalanobis) %>% 
mutate(MahalanobisRound = round(Mahalanobis,2))
# Calculating quantiles for each sample (for gene)
for( gene_ in colnames(horizontalTpm)[-1] ){
    # For each gene:
    # Filter gene associated information
    aux = hlaExp %>% filter(gene_name == gene_)
    cdf_ = ecdf(aux$tpm)
    aux[,gsub('-','',gene_)] = paste0(round(100*cdf_(aux$tpm),2))
    filterSamples = aux %>% filter(sampleid %in% samplesDf$sampleid) %>% select(sampleid,gsub('-','',gene_))
    samplesDf = merge(samplesDf ,filterSamples)

}

# Saving tables
samplesDf = samplesDf %>% arrange(-Mahalanobis)
outliersPt1 = samplesDf[,c(1,3:8)]
outliersPt2 = samplesDf[,c(1,3,9:12)]

write.table(outliersPt1,'/raid/genevol/users/lucas/heritability/00.Descriptive/Tables/multiOutPt1.txt',sep = '&',quote = F , row.names = F,col.names = T)
write.table(outliersPt2,'/raid/genevol/users/lucas/heritability/00.Descriptive/Tables/multiOutPt2.txt',sep = '&',quote = F , row.names = F,col.names = T)


# Comparing distributions of gene expressions with normal distribution
normalComparisonDf = as.data.frame(NULL)
for( idx_ in 1:length(genesLoop) ){

    # for each gene
    print(gene_)
    gene_ = genesLoop[idx_]
    # Filter gene associated information
    aux = hlaExp %>% filter(gene_name == gene_)
    # Calculate metrics
    ks_ = ks.test(aux$tpm,'pnorm',mean(aux$tpm),sd(aux$tpm))
    ad_ = nortest::ad.test(aux$tpm)
    check_ = length(unique(aux$tpm)) == 445
    normalComparisonDf[idx_,'Gene'] = gene_
    normalComparisonDf[idx_,'KS'] = round(100*ks_$p.value,2)
    normalComparisonDf[idx_,'AD'] = round(100*ad_$p.value,2)
    normalComparisonDf[idx_,'Ties'] = ifelse(check_ == T , 'No' , 'Yes')
    normalComparisonDf[idx_,'Skewness'] = moments::skewness(aux$tpm)
    
}

# Saving table
write.table(normalComparisonDf,paste0(rootSaveTables,'normalityTest.txt'), sep = '&', quote = F,row.names = F)

# Checking if data is discrete
for( idx_ in 1:length(genesLoop) ){
    
    print(gene_)
    gene_ = genesLoop[idx_]
    aux = hlaExp %>% filter(gene_name == gene_)
    print(length(unique(aux$tpm)))
    
}
# It is not


# Plot phenotypes by sex (boxplot)
geneExpressionSex = hlaExp %>% ggplot(aes(x = sex , y=tpm , fill = sex)) +
geom_boxplot(alpha = .8) +
facet_wrap(~gene_name , scale = 'free_y') +
theme_bw() +
theme( 
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
labs(x = 'Gene',y= "Transcriptions per million" , fill = 'Gene')

ggsave(paste0(rootSave,'densityPhenoSex',".png"), geneExpressionSex, bg = "transparent",width=10, height=8, dpi=300)

# Homoscedasticity by sex
bartlettSex = as.data.frame(NULL)
for( idx_ in 1:length(genesLoop) ){
    
    gene_ = genesLoop[idx_]
    aux = hlaExp %>% filter(gene_name == gene_)
    bart = bartlett.test(aux$tpm ~aux$sex)
    bartlettSex[idx_,'Gene'] = gene_
    bartlettSex[idx_,'pValue'] = bart$p.value
    aux2 = aux %>% filter(pop != 'YRI')
    bart2 = bartlett.test(aux2$tpm ~aux2$sex)
    bartlettSex[idx_,'pValueNA'] = bart2$p.value
    
}
confidence = .95
correctedAlpha = (1-confidence)/(2*length(genesLoop))
bartlettSex = bartlettSex %>% mutate(rejectH0 = ifelse(pValue < correctedAlpha , 'Yes' , 'No'),
                                      rejectH0NAfr = ifelse(pValueNA < correctedAlpha , 'Yes' , 'No')
                                    )

# Plot phenotypes by pop (boxplot)
geneExpressionPop = hlaExp %>% ggplot(aes(x = pop , y=tpm , fill = pop)) +
geom_boxplot(alpha = .8) +
facet_wrap(~gene_name , scale = 'free_y') +
theme_bw() +
theme( 
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
labs(x = 'Gene',y= "Transcriptions per million" , fill = 'Gene')

ggsave(paste0(rootSave,'densityPhenoPop',".png"), geneExpressionPop, bg = "transparent",width=10, height=8, dpi=300)

# Homoscedasticity by pop
bartlettPop = as.data.frame(NULL)
for( idx_ in 1:length(genesLoop) ){
    
    gene_ = genesLoop[idx_]
    aux = hlaExp %>% filter(gene_name == gene_)
    bart = bartlett.test(aux$tpm ~aux$pop)
    bartlettPop[idx_,'Gene'] = gene_
    bartlettPop[idx_,'pValue'] = bart$p.value
    
}
bartlettPop = bartlettPop %>% mutate(rejectH0 = ifelse(pValue < correctedAlpha , 'Yes' , 'No'))

# Plot phenotypes by sex and pop
geneExpressionTri = hlaExp %>% ggplot(aes(x = pop , y=tpm , fill = sex)) +
geom_boxplot(alpha = .8) +
facet_wrap(~gene_name , scale = 'free_y') +
theme_bw() +
theme( 
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
labs(x = 'Gene',y= "Transcriptions per million" , fill = 'Gene')

ggsave(paste0(rootSave,'densityPhenoTri',".png"), geneExpressionTri, bg = "transparent",width=10, height=8, dpi=300)

bartlettTri = as.data.frame(NULL)
for( idx_ in 1:length(genesLoop) ){
    
    gene_ = genesLoop[idx_]
    aux = hlaExp %>% filter(gene_name == gene_)
    aux$Join = paste0(aux$sex , aux$pop)
    bart = bartlett.test(aux$tpm ~aux$Join)
    bartlettTri[idx_,'Gene'] = gene_
    bartlettTri[idx_,'pValue'] = bart$p.value
    
}
bartlettTri = bartlettTri %>% mutate(rejectH0 = ifelse(pValue < correctedAlpha , 'Yes' , 'No'))
