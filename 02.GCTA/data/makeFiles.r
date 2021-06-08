listGeuvadis = read.table('/raid/genevol/heritability/samples.txt',header = F,stringsAsFactors = F)
colnames(listGeuvadis) = "samples"

hlaExp = readr::read_tsv("/raid/genevol/heritability/hla_expression.tsv") %>% rename(sampleid = subject_id) %>% filter(sampleid %in% listGeuvadis$samples)
for(gene_ in c('HLA-A', 'HLA-B', 'HLA-C', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DRA', 'HLA-DRB1')){
    aux_ = hlaExp  %>% filter(gene_name == gene_) %>% mutate(samp = sampleid) %>% select(samp,sampleid,tpm)
    write.table(aux_ , paste0('/raid/genevol/users/lucas/heritability/02.GCTA/data/gene' , gene_,'.txt'),col.names = F , row.names = F , quote = F)
}
cov1 = hlaExp  %>% filter(gene_name == 'HLA-A') %>% mutate(samp = sampleid) %>% select(samp,sampleid,sex,lab)
write.table(cov1,'/raid/genevol/users/lucas/heritability/02.GCTA/data/envir1.txt',col.names = F , row.names = F , quote = F)
cov2 = hlaExp  %>% filter(gene_name == 'HLA-A') %>% mutate(samp = sampleid) %>% select(samp,sampleid,lab)
write.table(cov2,'/raid/genevol/users/lucas/heritability/02.GCTA/data/envir2.txt',col.names = F , row.names = F , quote = F)

