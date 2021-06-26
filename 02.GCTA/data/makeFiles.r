listGeuvadis = read.table('/raid/genevol/heritability/samples.txt',header = F,stringsAsFactors = F)
colnames(listGeuvadis) = "samples"

hlaExp = readr::read_tsv("/raid/genevol/heritability/hla_expression.tsv") %>% rename(sampleid = subject_id) %>% filter(sampleid %in% listGeuvadis$samples)
hlaExpNAfr = hlaExp %>% filter(pop != 'YRI')
for(gene_ in c('HLA-A', 'HLA-B', 'HLA-C', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DRA', 'HLA-DRB1')){
    aux_ = hlaExp  %>% filter(gene_name == gene_) %>% mutate(samp = sampleid) %>% select(samp,sampleid,tpm)
    write.table(aux_ , paste0('/raid/genevol/users/lucas/heritability/02.GCTA/data/gene' , gene_,'.txt'),col.names = F , row.names = F , quote = F)
    aux2_ = hlaExpNAfr  %>% filter(gene_name == gene_) %>% mutate(samp = sampleid) %>% select(samp,sampleid,tpm)
    write.table(aux2_ , paste0('/raid/genevol/users/lucas/heritability/02.GCTA/data/geneNAfr' , gene_,'.txt'),col.names = F , row.names = F , quote = F)

    for(version_ in c(1,2,3)){
        dir.create(paste0('/raid/genevol/users/lucas/heritability/03.Heritability/Results/V',version_,'/',gene_))
    }
}

cov0 = hlaExp  %>% filter(gene_name == 'HLA-A') %>% mutate(samp = sampleid) %>% select(samp,sampleid,sex,lab,pop)
write.table(cov0,'/raid/genevol/users/lucas/heritability/02.GCTA/data/envir0.txt',col.names = F , row.names = F , quote = F)
cov1 = hlaExp  %>% filter(gene_name == 'HLA-A') %>% mutate(samp = sampleid) %>% select(samp,sampleid,sex,lab)
write.table(cov1,'/raid/genevol/users/lucas/heritability/02.GCTA/data/envir1.txt',col.names = F , row.names = F , quote = F)
cov2 = hlaExp  %>% filter(gene_name == 'HLA-A') %>% mutate(samp = sampleid) %>% select(samp,sampleid,lab)
write.table(cov2,'/raid/genevol/users/lucas/heritability/02.GCTA/data/envir2.txt',col.names = F , row.names = F , quote = F)

