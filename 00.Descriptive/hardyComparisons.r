suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(gtable)))

manhattanPltDf = function(df_,samp_){
#adapted from:https://www.r-graph-gallery.com/101_Manhattan_plot.html
    dfPrep <- df_ %>% 
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=as.numeric(max(BP))) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(df_, ., by=c("CHR"="CHR")) %>%
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate(BPcum=BP+tot , Samp = samp_)

  axisdf = dfPrep %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )%>%
  ungroup() %>%
  mutate(Samp = samp_)

  return(list(dfPrep,axisdf))

}

manhattanPlt = function(dfPrep_,axisdf_,cutPoint_,type_){
#adapted from:https://www.r-graph-gallery.com/101_Manhattan_plot.html

  if(type_ == 'hwe'){
    plt0 = ggplot(dfPrep_, aes(x=BPcum, y=-log10(P)))
  } else{
    plt0 = ggplot(dfPrep_, aes(x=BPcum, y=MAF))
  }
    plt = plt0 +
    # Show all points
    geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf_$CHR, breaks= axisdf_$center ) +
    scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis

    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
    geom_hline(yintercept = cutPoint_ , col = 'red') +
    facet_wrap(~Samp , ncol = 1)
    
    return(plt)
}

savePlt = function(df_,sample_,root_,name_){

    dfPlot = df_ %>% filter(Sample == sample_)
    p = manhattanPlt(dfPlot,1e-4)
    ggsave(paste0(root_,name_,".png"), p, bg = "transparent")
    return(p)

}

root = '/scratch/genevol/users/lucas/'
rootSave = '/raid/genevol/users/lucas/heritability/00.Descriptive/Images/Hardy/'
hweDf = NULL
mafDf = NULL
refDf = NULL

# Reading HWE files and SNP's mappings

for(i in 1:22){

    print(paste0('Chromosome: ',i))
    # Reading p-value from hwe test for each SNP from i-th chromosome
    # Hwe test from individuals whose phenotype is known, considering only europeans (n=358)
    auxNafrHwe = data.table::fread(paste0(root,'Hardy_NAfr_chr',i,'.hwe'), header = T,select = c('CHR' , 'SNP' , 'P')) %>% 
    mutate(Sample = 'NAfr') # Labelling df
    # Maf from individuals whose phenotype is known, considering only europeans (n=358)
    auxNafrMaf = data.table::fread(paste0(root,'Maf/Hardy_NAfr_chr',i,'.frq'), header = T,select = c('CHR' , 'SNP' , 'MAF')) %>% 
    mutate(Sample = 'NAfr') # Labelling df
    # Hwe test from individuals whose phenotype is known (n=445)
    auxGeuvadisHwe = data.table::fread(paste0(root,'Hardy_Geuvadis_chr',i,'.hwe'), header = T,select = c('CHR' , 'SNP' , 'P')) %>% 
    mutate(Sample = 'Geuvadis') # Labelling df
    # Maf from individuals whose phenotype is known (n=445)
    auxGeuvadisMaf = data.table::fread(paste0(root,'Maf/Hardy_Geuvadis_chr',i,'.frq'), header = T,select = c('CHR' , 'SNP' , 'MAF')) %>% 
    mutate(Sample = 'Geuvadis') # Labelling df
    # Hwe test from all genotyped individuals (n=2504)
    auxFullHwe = data.table::fread(paste0(root,'Hardy_Full_chr',i,'.hwe'), header = T,select = c('CHR' , 'SNP' , 'P')) %>% 
    mutate(Sample = 'Full') # Labelling df
    # Maf from all genotyped individuals (n=2504)
    auxFullMaf = data.table::fread(paste0(root,'Maf/Hardy_Full_chr',i,'.frq'), header = T,select = c('CHR' , 'SNP' , 'MAF')) %>% 
    mutate(Sample = 'Full') # Labelling df
    # Reading SNP's position
    auxRef = data.table::fread(paste0(root,'basicInfoChr',i,'.txt'), header = F)
    colnames(auxRef) = c('CHR','SNP','BP') 
    
    # Appending tests
    hweDf = rbind(hweDf , rbind(auxNafrHwe,rbind(auxGeuvadisHwe,auxFullHwe)))
    # Appending MAF
    mafDf = rbind(mafDf , rbind(auxNafrMaf,rbind(auxGeuvadisMaf,auxFullMaf)))
    # Appending mapping from SNPs
    refDf = rbind(refDf , auxRef)

}

# merge between hwe test and SNP positions
finalPltDfHwe = merge(hweDf,refDf)
# merge between MAF and SNP positions
finalPltDfMaf = merge(mafDf,refDf)


dfPltManhattanHwe = NULL
dfPltAxisManhattanHwe = NULL
dfPltManhattanMaf = NULL
dfPltAxisManhattanMaf = NULL

filters_ = c('NAfr','Geuvadis','Full')
names_ = c('Only europeans' , 'Known phenotypes' , 'All genotyped individuals')
for(idx_ in 1:3){

    print(paste0('HWE, sample:',names_[idx_]))
    auxFiltHwe = finalPltDfHwe %>% filter(Sample == filters_[idx_])
    auxListHwe = manhattanPltDf(auxFiltHwe,names_[idx_])
    dfPltManhattanHwe = rbind(dfPltManhattanHwe , auxListHwe[[1]])
    dfPltAxisManhattanHwe = rbind(dfPltAxisManhattanHwe , auxListHwe[[2]])
    auxListHwe = NULL
    auxFiltHwe = NULL

    print(paste0('MAF, sample:',names_[idx_]))
    auxFiltMaf = finalPltDfMaf %>% filter(Sample == filters_[idx_])
    auxListMaf = manhattanPltDf(auxFiltMaf,names_[idx_])
    dfPltManhattanMaf = rbind(dfPltManhattanMaf , auxListMaf[[1]])
    dfPltAxisManhattanMaf = rbind(dfPltAxisManhattanMaf , auxListMaf[[2]])
    auxListMaf = NULL

}

Sys.time()
plotMHwe = manhattanPlt(dfPltManhattanHwe,dfPltAxisManhattanHwe,4,'hwe')
ggsave(paste0(rootSave,'allFramesHwe',".png"), plotMHwe, bg = "transparent",width=12, height=6, dpi=300)
Sys.time()

dfHwe = NULL
cuts = c(.05 , 1e-2 , 1e-3 , 1e-4 , 1e-5 , 1e-6 , 1e-7 , 1e-8 )
totSnps = dfPltManhattanMaf %>% group_by(Sample,CHR) %>% summarise(nTot=n()) %>% ungroup()
for(cut_ in cuts){
  print(cut_)
  dfCuts = merge(
    finalPltDfHwe %>% filter(P <= cut_) %>% group_by(Sample,CHR) %>% summarise(n=n()) %>% ungroup() %>% mutate(cut = cut_)
    ,
    totSnps
  )
  dfHwe = rbind(dfHwe,dfCuts)
}


dfHwe = dfHwe %>% mutate(freq = n/nTot)

dfMaf = NULL
cuts = c(.001 , .005 , (1:5)/100)

for(cut_ in cuts){
  print(cut_)
  dfCuts = merge(
    dfPltManhattanMaf %>% filter(MAF <= cut_) %>% group_by(Sample,CHR) %>% summarise(n=n()) %>% ungroup() %>% mutate(cut = cut_)
    ,
    totSnps
  )
  dfMaf = rbind(dfMaf,dfCuts)
}

dfMaf = dfMaf %>% mutate(freq = n/nTot)

knitr::kable(dfMaf)

plotMafLog = dfMaf %>% ggplot(aes(x = CHR,y=log10(n) , colour = as.factor(cut))) + 
geom_point() +
geom_line() + 
facet_wrap(~Sample,ncol = 1) +
theme_bw() +
theme( 
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
labs(x = 'Chromosome',y= expression(Frequency ~ (log[10]) ) , colour = 'Cutpoint')

plotMafRelativeFr = dfMaf %>% ggplot(aes(x = CHR,y = freq , colour = as.factor(cut))) + 
geom_point() +
geom_line() + 
facet_wrap(~Sample,ncol = 1) +
theme_bw() +
theme( 
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
labs(x = 'Chromosome',y= "Relative frequency" , colour = 'Cutpoint')

grobLog = ggplotGrob(plotMafLog)
grobLog = gtable_add_cols(grobLog,unit(0,'mm'))
grobRel = ggplotGrob(plotMafRelativeFr)
grobRel = gtable_add_cols(grobRel,unit(0,'mm'))
bindGrobs = rbind(grobLog,grobRel , size = 'first')


ggsave(paste0(rootSave,'allFramesMafx',".png"), bindGrobs, bg = "transparent",width=10, height=8, dpi=300)

