library(coxme)
library(dplyr)
library(yaml)



##########################################################################################
# Functions
##########################################################################################
matrix_prep_ = function(matrix_){
    matrix_ = matrix_ %>% select(-V1)
    colnames(matrix_) = subject_ids$subject_id
    rownames(matrix_) = subject_ids$subject_id
    return(matrix_)
}

k_fold_samples_ = function(df_, k_){
    set.seed(75952)
    total_sample_ = nrow(df_)
    n_each_fold_ = floor(total_sample_/10)
    idx_samples_ = list()
    df_$set_ = 0
    for(i_ in c(1:k_)){
        if (i_ < k_){
            aux_ = df_ %>% 
                    filter(set_ == 0)
            sample_ = aux_ %>% filter(subject_id %in% sample(aux_$subject_id, n_each_fold_))
        }
        df_[df_$subject_id %in% sample_$subject_id,'set_'] = i_
    }
    folds_ = list()
    for(i_ in c(1:k_)){
        folds_[[i_]] = df_[df_$set_ != i_,'subject_id']
    }
    return(folds_)
}

fit_model_ = function(df_, varlist_){
    new_varlist_ = list()
    for (set_ in names(varlist_)){
        matrix_ = varlist_[[set_]]
        matrix_ = matrix_[rownames(matrix_) %in% df_$subject_id, colnames(matrix_) %in% df_$subject_id]
        new_varlist_[[set_]] = as.matrix(matrix_)
    }
    fit1 = lmekin(
                as.formula(paste("tpm~",formula_)), 
                data = df_, 
                random=~1|df_$subject_id, 
                varlist=new_varlist_, 
                vinit=2, 
                method = "REML"
            )
    sig2e1 = fit1$sigma^2 ## Captura a variância do erro
    sig2g1 = as.numeric(fit1$vcoef) ## Captura a variância genética (poligênica)
    h1 = sig2g1/(sig2g1+sig2e1) # Cálculo da herdabilidade
    return(list(h1 = h1, sig2e1 = sig2e1, sig2g1 = sig2g1))
}



estimate_h2 = function(list_){
    dataset_ = list_[['dataset_']]
    df_ = list_[['df_']]
    varlist_ = list_[['varlist_']]
    genes_= list_[['genes_']]
    formula_ = list_[['formula_']]
    return_ = list()
    for(gene_ in genes_){
        vector_h1_ = c()
        phen_data = dataset_ %>% 
            filter(gene_name == gene_) %>% 
            select(tpm, subject_id)
        df_complete_fit_ = merge(df_, phen_data, by = c('subject_id'))
        # folds_fit_ = k_fold_samples_(df_complete_fit_, 10)
        B = 10
        frac_ = .9
        vector_h1_ = c()
        set.seed(75952)
        for(fold_ in c(1:(B))){
            # df_sample_ = df_complete_fit_ %>% filter(subject_id %in% folds_fit_[[fold_]])
            df_sample_ = df_complete_fit_ %>% filter(subject_id %in% sample(df_complete_fit_$subject_id, floor(nrow(df_complete_fit_) * frac_)))
            vector_h1_ = c(vector_h1_, fit_model_(df_sample_, varlist_)[['h1']])
        }
        print(vector_h1_)
        df_sample_ = df_complete_fit_ %>% filter(subject_id %in% sample(df_complete_fit_$subject_id, floor(nrow(df_complete_fit_) * frac_)))
        fit_ = fit_model_(df_sample_, varlist_)
        h1 = fit_[['h1']]
        sig2e1 = fit_[['sig2e1']]
        sig2g1 = fit_[['sig2g1']]
        ELRT = -sum(log(1+h1*(eigen(varlist_[[1]])$values-1)))
        p_value_ = (1 - pchisq(ELRT, 1)) / 2
        sd_h1 = sd(vector_h1_)
        return_ = list(
                        sig2e1=sig2e1,
                        sig2g1=sig2g1,
                        h1=h1,
                        ELRT = ELRT,
                        p_value_=p_value_, 
                        sd_h1 = sd_h1
                        )
    }
    return(return_)
}



##########################################################################################
# Process
##########################################################################################
path_tmp_ = commandArgs(TRUE)[1]
formula_ = commandArgs(TRUE)[2]
set_ = commandArgs(TRUE)[3]
matrix_name_ = commandArgs(TRUE)[4]
gene_ = commandArgs(TRUE)[5]
yml_ = read_yaml(paste0('conf/',tail(strsplit(path_tmp_, '/')[[1]],1), '/parameters.yml'))


print(path_tmp_)
print(formula_)
print(set_)
print(matrix_name_)
print(gene_)



dataset_ = read.table("data/01_raw/hla_expression.tsv",sep = '\t', header = T)
# All fam files are equal - chose chromosome 22 arbitrarily
subject_ids = read.table(paste0(path_tmp_,"/pruned_chr22.fam"), sep = ' ') %>% 
                    select(V1) %>% 
                    rename(subject_id = V1)

# The dataset is stacked
# In order to read sex, lab and pop, it is only needed to filter a gene
# (in this case was HLA-A) and then drop duplicates
fixed_covs = dataset_ %>% 
            filter(gene_name == 'HLA-A') %>% 
            unique() %>% 
            select(-c(gene_name, tpm))

df_fit = merge(subject_ids, fixed_covs, all.x = T, all.y = F, by = 'subject_id')
varlist_ = list()
for(matrix_ in yml_$params_run$matrices$sets_[set_]){
    matrix_df_ = read.table(paste0(path_tmp_,'/',matrix_name_,'_',set_,'_',names(matrix_),'_correction_non_negative.txt'), sep = '|')
    varlist_[[names(matrix_)]] = matrix_prep_(matrix_df_)
}


list_ = list(
                df_ = df_fit
                , dataset_ = dataset_
                , varlist_ = varlist_
                , genes_ = gene_
                , formula_ = formula_
            )
result_ = as.data.frame(estimate_h2(list_))
write.table(result_,
            paste0(
                    path_tmp_,'/',
                    'resultfit_',
                    formula_,'_',
                    set_, '_',
                    matrix_name_,'_',
                    gene_,'.txt'
                    )
            )

