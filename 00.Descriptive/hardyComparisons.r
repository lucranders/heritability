suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(ggplot2)))

root = '/scratch/genevol/users/lucas/'

hweDf = NULL
refDf = NULL

for(i in 1:22){

    print(paste0('Chromosome: ',i))
    # Reading p-value from hwe test for each SNP from i-th chromosome
    # Hwe test from individuals whose phenotype is known, considering only europeans (n=358)
    auxNafr = data.table::fread(paste0(root,'Hardy_NAfr_chr',i,'.hwe'), header = T,select = c('CHR' , 'SNP' , 'P')) %>% 
    mutate(Sample = 'NAfr') # Labelling df
    # Hwe test from individuals whose phenotype is known (n=445)
    auxGeuvadis = data.table::fread(paste0(root,'Hardy_Geuvadis_chr',i,'.hwe'), header = T,select = c('CHR' , 'SNP' , 'P')) %>% 
    mutate(Sample = 'Geuvadis') # Labelling df
    # Hwe test from all genotyped individuals (n=2504)
    auxFull = data.table::fread(paste0(root,'Hardy_Full_chr',i,'.hwe'), header = T,select = c('CHR' , 'SNP' , 'P')) %>% 
    mutate(Sample = 'Full') # Labelling df
    # Reading SNP's position
    auxRef = data.table::fread(paste0(root,'basicInfoChr',i,'.txt'), header = F)
    colnames(auxRef) = c('CHR','SNP','BP') 
    
    # Appending tests
    hweDf = rbind(hweDf , rbind(auxNafr,rbind(auxGeuvadis,auxFull)))
    # Appending mapping from SNPs
    refDf = rbind(refDf , auxRef)

}

# merge between hwe test and SNP positions
finalPltDf = merge(hweDf,refDf)

qqman::manhattan(finalPltDf %>% filter(Sample == 'NAfr'),suggestiveline =  -log10(1e-04),main = "HW disequilibrium, considering only europeans which phenotype is known (n=358)")
p = recordPlot()
g = grid::grid.grabExpr(grid::grid.echo(p))
ggsave("/raid/genevol/users/lucas/heritability/00.Descriptive/Images/Hardy/nAfr.png", g, bg = "transparent")

qqman::manhattan(finalPltDf %>% filter(Sample == 'Geuvadis'),suggestiveline =  -log10(1e-04), main = "HW disequilibrium, considering only individuals which phenotype is known (n=445)")
p = recordPlot()
g = grid::grid.grabExpr(grid::grid.echo(p))
ggsave("/raid/genevol/users/lucas/heritability/00.Descriptive/Images/Hardy/geuvadis.png", g, bg = "transparent")

qqman::manhattan(finalPltDf %>% filter(Sample == 'Full'),suggestiveline =  -log10(1e-04), main = "HW disequilibrium, all genotyped individuals (n=2504)")
p = recordPlot()
g = grid::grid.grabExpr(grid::grid.echo(p))
ggsave("/raid/genevol/users/lucas/heritability/00.Descriptive/Images/Hardy/Full.png", g, bg = "transparent")

table(finalPltDf$CHR)

