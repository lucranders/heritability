suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(ggplot2)))

rootSave = '/raid/genevol/users/lucas/heritability/00.Descriptive/Images/Phenotype/'
rootSaveTables = '/raid/genevol/users/lucas/heritability/00.Descriptive/Tables/'
listGeuvadis = read.table('/raid/genevol/heritability/samples.txt',header = F)
colnames(listGeuvadis) = "samples"
hlaExp = readr::read_tsv("/raid/genevol/heritability/hla_expression.tsv") %>% rename(sampleid = subject_id) %>% filter(sampleid %in% listGeuvadis$samples)

##################################################################################################################################
# Distribution of genes and symmetry analysis
##################################################################################################################################


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


##################################################################################################################################
# Multivariate analysis - checking for possible outliers - Mahalanobis
##################################################################################################################################

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

dfMahalanobis = function(horizontalExpDf,prob_=.993,df_=9){

    # Filter columns of interest - gene expressions and samples
    horizontalTpm = horizontalExpDf[,grepl('HLA|sample',colnames(horizontalExpDf))] 

    # Calculate means and covariance matrix
    meanExp = colMeans(horizontalTpm[,2:10])
    covExp = cov(horizontalTpm[,2:10])
    # Calculate Mahalanobis distance
    mahalanobisD = mahalanobis(horizontalTpm[,2:10],center = meanExp,cov = covExp)
    horizontalExpDf$Mahalanobis = mahalanobisD
    horizontalExpDf$idx = 1:nrow(horizontalExpDf)
    # Defining threshold on 99.9% chi-squared (9 df) quantile 
    cutMahalanobis = qchisq(prob_, df = df_)
    # Create dummy variable to define colours in the plot
    horizontalExpDf = horizontalExpDf %>% mutate(threshold = ifelse(Mahalanobis > cutMahalanobis , 'yes' , 'no'))

    return(horizontalExpDf)

}

plotMahalanobis = function(horizontalExpDf,prob_=.993,df_=9){

    cutMahalanobis = qchisq(prob_, df = df_)
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

    return(plotMahalanobis)

}
# ggsave(paste0(rootSave,'mahalanobis',".png"), plotMahalanobis, bg = "transparent",width=10, height=8, dpi=300)

plotMahalanobisQQ = function(horizontalExpDf,prob_=.993,df_=9){

    cutMahalanobis = qchisq(prob_, df = df_)
    # Sorting Mahalanobis values
    dfOrderedMahalanobis = as.data.frame( sort ( horizontalExpDf$Mahalanobis) )
    # Defining quantile values to compare with Chi-squared distribution
    dfOrderedMahalanobis$idx = ((1:nrow(dfOrderedMahalanobis)) - .5)/nrow(dfOrderedMahalanobis)
    dfOrderedMahalanobis$chisqQ = qchisq(dfOrderedMahalanobis$idx, df = df_)
    # Selecting only observed mahalanobis and expected chi-squared value
    dfOrderedMahalanobisPlt = dfOrderedMahalanobis[,c(1,3)]
    colnames(dfOrderedMahalanobisPlt) = c('Observed','Expected')
    # Label high values
    dfOrderedMahalanobisPlt = dfOrderedMahalanobisPlt %>% mutate(Threshold = ifelse(Observed > cutMahalanobis , 'yes' , 'no'))

    # Creating confidence bands
    intervalFull = confidenceBandsMahalanobis(n = nrow(dfOrderedMahalanobisPlt) , idxInf = 50 , idxSup = 9950 , qtTrials = 10000)
    dfOrderedMahalanobisPlt$LimInf = intervalFull[[1]]
    dfOrderedMahalanobisPlt$LimSup = intervalFull[[2]]

    # Quantile-quantile plot - comparing Mahalanobis distance with Chi-squared distribution
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

    return(plotQQFull)
    
}

# Transforming original dataset in horizontal position to analyse
dfMahalanobisAnalysis = tidyr::spread(hlaExp,gene_name,tpm)
dfMahalanobisEurAnalysis = tidyr::spread(hlaExp,gene_name,tpm) %>% filter(pop != 'YRI')

# Mahalanobis with all available data
dfMahalanobisFull = dfMahalanobis(dfMahalanobisAnalysis)
# Mahalanobis, excluding possible outliers
dfMahalanobisOut = dfMahalanobis(dfMahalanobisFull %>% filter(threshold == 'no'))
# Mahalanobis with only europeans
dfMahalanobisEur = dfMahalanobis(dfMahalanobisEurAnalysis)
# Mahalanobis with only europeans, excluding possible outliers
dfMahalanobisEurOut = dfMahalanobis(dfMahalanobisEur %>% filter(threshold == 'no'))

p11 = plotMahalanobis(dfMahalanobisFull)
p12 = plotMahalanobis(dfMahalanobisEur)
p21 = plotMahalanobis(dfMahalanobisOut)
p22 = plotMahalanobis(dfMahalanobisEurOut)

pltMahalanobis = ggpubr::ggarrange( 
                    p11, p12, p21 , p22 ,
                    labels = c("Entire", "Europeans",'Entire - Filtered','Europeans - Filtered'),
                    ncol = 2, nrow = 2
                    )

p11 = plotMahalanobisQQ(dfMahalanobisFull)
p12 = plotMahalanobisQQ(dfMahalanobisEur)
p21 = plotMahalanobisQQ(dfMahalanobisOut)
p22 = plotMahalanobisQQ(dfMahalanobisEurOut)

pltMahalanobis = ggpubr::ggarrange( 
                    p11, p12, p21 , p22 ,
                    labels = c("Entire", "Europeans",'Entire - Filtered','Europeans - Filtered'),
                    ncol = 2, nrow = 2
                    )


ggsave(paste0(rootSave,'mahalanobisCompilationQQ',".png"), pltMahalanobis, bg = "transparent",width=10, height=8, dpi=300)

# Selecting possible outliers
samplesDf = dfMahalanobisFull %>% filter(threshold == 'yes') %>% select(sampleid,pop,lab,Mahalanobis) %>% 
mutate(MahalanobisRound = round(Mahalanobis,2))

samplesDfEur = dfMahalanobisEurOut %>% select(sampleid,Mahalanobis) %>% 
mutate(MahalanobisRound = round(Mahalanobis,2))

# Calculating quantiles for each sample (for gene)
for( gene_ in colnames(dfMahalanobisOut)[grepl('HLA',colnames(dfMahalanobisOut))] ){
    # For each gene:
    # Filter gene associated information
    print(gene_)
    dfOrderedTpm = as.data.frame( dfMahalanobisAnalysis %>% select(gene_) %>% arrange(get(gene_)) )
    dfOrderedTpmEur = as.data.frame( dfMahalanobisEurAnalysis %>% select(gene_) %>% arrange(get(gene_)) )
    # Defining quantile values to compare with Chi-squared distribution
    dfOrderedTpm$idx = ((1:nrow(dfOrderedTpm)) - .5)/nrow(dfOrderedTpm)
    dfOrderedTpmEur$idx = ((1:nrow(dfOrderedTpmEur)) - .5)/nrow(dfOrderedTpmEur)
    colnames(dfOrderedTpm) = c(gene_,'idx')
    colnames(dfOrderedTpmEur) = c(gene_,'idx')
    cutPoint = dfOrderedTpm  %>% filter(idx > .993) %>% summarise(cut_ = min(get(gene_)))
    cutPointEur = dfOrderedTpmEur  %>% filter(idx > .993) %>% summarise(cut_ = min(get(gene_)))

    aux = merge(dfMahalanobisAnalysis[,c('sampleid',gene_)],dfOrderedTpm)
    auxEur = merge(dfMahalanobisEurAnalysis[,c('sampleid',gene_)],dfOrderedTpmEur)
    aux[,gsub('-','',gene_)] = ifelse(aux[,gene_] >= cutPoint$cut_ , paste0('cellcolor{col2}', aux$idx) , aux$idx)
    auxEur[,gsub('-','',gene_)] = ifelse(auxEur[,gene_] >= cutPointEur$cut_ , paste0('cellcolor{col2}', auxEur$idx) , auxEur$idx)
    filterSamples = aux %>% filter(sampleid %in% samplesDf$sampleid) %>% select(sampleid,gsub('-','',gene_))
    filterSamplesEur = auxEur %>% filter(sampleid %in% samplesDfEur$sampleid) %>% select(sampleid,gsub('-','',gene_))
    samplesDf = merge(samplesDf ,filterSamples)
    samplesDfEur = merge(samplesDfEur ,filterSamplesEur)

}

# Saving tables
samplesDf = samplesDf %>% arrange(-Mahalanobis)
samplesDfEur = samplesDfEur %>% arrange(-Mahalanobis)
for (idx in 1:nrow(samplesDf)){
    samplesDf[idx,'count'] = sum(grepl('cellcolor',samplesDf[idx,]))
}

for (idx in 1:nrow(samplesDfEur)){
    samplesDfEur[idx,'count'] = sum(grepl('cellcolor',samplesDfEur[idx,]))
}

outliersPt1 = samplesDf[,c(1,3:8)]
outliersPt2 = samplesDf[,c(1,3,9:13)]

outliersPt1Eur = samplesDfEur[,c(1,3:8)]
outliersPt2Eur = samplesDfEur[,c(1,3,9:13)]


write.table(outliersPt1,'/raid/genevol/users/lucas/heritability/00.Descriptive/Tables/multiOutPt1.txt',sep = '&',quote = F , row.names = F,col.names = T)
write.table(outliersPt2,'/raid/genevol/users/lucas/heritability/00.Descriptive/Tables/multiOutPt2.txt',sep = '&',quote = F , row.names = F,col.names = T)

write.table(outliersPt1Eur,'/raid/genevol/users/lucas/heritability/00.Descriptive/Tables/multiOutEurPt1.txt',sep = '&',quote = F , row.names = F,col.names = T)
write.table(outliersPt2Eur,'/raid/genevol/users/lucas/heritability/00.Descriptive/Tables/multiOutEurPt2.txt',sep = '&',quote = F , row.names = F,col.names = T)


##################################################################################################################################
# Multivariate analysis - biplots
##################################################################################################################################


# Standardize gene expressions
biplotPlt = function(covars_,sampleId_,cor_,dfThreshold_,scaleBip_) {

    biplotFit = princomp(x = covars_, cor = cor_)
    biplotDf = data.frame(sampleid = sampleId_ , biplotFit$scores)
    # Merge df with dummy variable based on Mahalanobis distance and threshold
    biplotDf = merge(biplotDf,dfThreshold_)
    pcLoads = data.frame(biplotFit$loadings[,1:9])
    
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
    ggrepel::geom_text_repel(data = pcLoads , aes(x=scaleBip*Comp.1, y=scaleBip*Comp.2, label=row.names(pcLoads) , size = 5, vjust=1) ) +
    # geom_text(data = pcLoads , aes(x=scaleBip*Comp.1, y=scaleBip*Comp.2, label=row.names(pcLoads) , size = 5, vjust=1) ) +
    labs(x = 'Principal Component 1' , y = 'Principal Component 2')

    return( plotBiplot )

}
# horizontalTpm[,2:10] = scale(horizontalTpm[,2:10],center = T , scale = T)
# Create biplot - PCA for XX'and X'X

covarsFull_ = dfMahalanobisAnalysis[,grepl('HLA',colnames(dfMahalanobisAnalysis))]
sampleFull_ = dfMahalanobisAnalysis[,grepl('sample',colnames(dfMahalanobisAnalysis))] 
pltFull = biplotPlt( covars_ = covarsFull_ , sampleId_ = sampleFull_ , cor_ = FALSE , dfThreshold_ = horizontalExpDf[,c('sampleid','threshold')] , scaleBip_ = 2500)

pltFull = biplotPlt( covars_ = horizontalTpm[,2:10] , sampleId_ = horizontalTpm[,1] , cor_ = FALSE , dfThreshold_ = horizontalExpDf[,c('sampleid','threshold')] , scaleBip_ = 2500)



biplotFitFilt = princomp(x = newMaha[,2:10], cor=FALSE)
biplotDfFilt = data.frame(sampleid = newMaha[,1] , biplotFitFilt$scores)
# Merge df with dummy variable based on Mahalanobis distance and threshold
biplotDfFilt = merge(biplotDfFilt,filterOut[,c('sampleid','threshold')])
pcLoadsFilt = data.frame(biplotFitFilt$loadings[,1:9])


# Plot biplot
plotBiplotFilt = biplotDfFilt %>% ggplot() + 
geom_text(data=subset(biplotDfFilt, threshold == 'yes'),
            aes(Comp.1,Comp.2,label=sampleid,colour = threshold),alpha = .7) +
geom_text(data=subset(biplotDfFilt, threshold == 'no'),
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
geom_segment(data = pcLoadsFilt , aes(x=0, y=0,xend = scaleBip*Comp.1 , yend = scaleBip*Comp.2), arrow=arrow(length=unit(0.1,"cm")), alpha=0.75) +
ggrepel::geom_text_repel(data = pcLoadsFilt , aes(x=scaleBip*Comp.1, y=scaleBip*Comp.2, label=row.names(pcLoads) , size = 5, vjust=1) ) +
# geom_text(data = pcLoadsFilt , aes(x=scaleBip*Comp.1, y=scaleBip*Comp.2, label=row.names(pcLoads) , size = 5, vjust=1) ) +
labs(x = 'Principal Component 1' , y = 'Principal Component 2')


# ggsave(paste0(rootSave,'biplotGeneExp',".png"), plotBiplot, bg = "transparent",width=10, height=8, dpi=300)
biplotPlts = ggpubr::ggarrange( 
                    plotBiplot, plotBiplotFilt,
                    labels = c("Full", "Filter"),
                    nrow = 2 , ncol = 1
                    )
ggsave(paste0(rootSave,'biplotsCompNAfr',".png"), biplotPlts, bg = "transparent",width=10, height=8, dpi=300)

##################################################################################################################################
# Univariate analysis - checking for heteroscedasticity and normality
##################################################################################################################################

outL = horizontalExpDf %>% filter(threshold == 'no') %>% select(sampleid)
# Comparing distributions of gene expressions with normal distribution
normalComparisonDf = as.data.frame(NULL)
for( idx_ in 1:length(genesLoop) ){

    # for each gene
    print(gene_)
    gene_ = genesLoop[idx_]
    # Filter gene associated information
    aux = hlaExp %>% filter(gene_name == gene_)
    auxOut = hlaExp %>% filter(gene_name == gene_ , sampleid %in% outL$sampleid)
    # Calculate metrics
    ks_ = ks.test(aux$tpm,'pnorm',mean(aux$tpm),sd(aux$tpm))
    ksOut_ = ks.test(auxOut$tpm,'pnorm',mean(auxOut$tpm),sd(auxOut$tpm))
    ad_ = nortest::ad.test(aux$tpm)
    adOut_ = nortest::ad.test(auxOut$tpm)
    check_ = length(unique(aux$tpm)) == nrow(aux)
    checkOut_ = length(unique(auxOut$tpm)) == nrow(auxOut)
    normalComparisonDf[idx_,'Gene'] = gene_
    normalComparisonDf[idx_,'KS'] = round(100*ks_$p.value,2)
    normalComparisonDf[idx_,'KSOut'] = round(100*ksOut_$p.value,2)
    normalComparisonDf[idx_,'AD'] = round(100*ad_$p.value,2)
    normalComparisonDf[idx_,'ADOut'] = round(100*adOut_$p.value,2)
    normalComparisonDf[idx_,'Ties'] = ifelse(check_ == T , 'No' , 'Yes')
    normalComparisonDf[idx_,'TiesOut'] = ifelse(checkOut_ == T , 'No' , 'Yes')
    normalComparisonDf[idx_,'Skewness'] = moments::skewness(aux$tpm)
    normalComparisonDf[idx_,'SkewnessOut'] = moments::skewness(auxOut$tpm)
    
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
