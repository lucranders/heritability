suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(ggplot2)))

rootSave = '/raid/genevol/users/lucas/heritability/00.Descriptive/Images/Phenotype/'
listGeuvadis = read.table('/raid/genevol/heritability/samples.txt',header = F)
colnames(listGeuvadis) = "samples"
hlaExp = readr::read_tsv("/raid/genevol/heritability/hla_expression.tsv") %>% rename(sampleid = subject_id) %>% filter(sampleid %in% listGeuvadis$samples)


# Plot phenotypes by gene (boxplot)
densityPlot = hlaExp %>% ggplot(aes(x = gene_name , y=tpm , fill = gene_name)) +
geom_boxplot(alpha = .8) +
theme_bw() +
theme( 
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
labs(x = 'Gene',y= "Transcriptions per million" , fill = 'Gene')

ggsave(paste0(rootSave,'densityPheno',".png"), densityPlot, bg = "transparent",width=10, height=8, dpi=300)

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
genesLoop = sort(unique(hlaExp$gene_name))
for( idx_ in 1:length(genesLoop) ){
    
    gene_ = genesLoop[idx_]
    aux = hlaExp %>% filter(gene_name == gene_)
    bart = bartlett.test(aux$tpm ~aux$sex)
    bartlettSex[idx_,'Gene'] = gene_
    bartlettSex[idx_,'pValue'] = bart$p.value
    
}
confidence = .95
correctedAlpha = (1-confidence)/(2*length(genesLoop))
bartlettSex = bartlettSex %>% mutate(rejectH0 = ifelse(pValue < correctedAlpha , 'Yes' , 'No'))

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
