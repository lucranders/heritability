from kedro.config import ConfigLoader
from plotnine import *
import numpy as np
import pandas as pd
import pickle

def matrix_to_df(params, id_data):
    matrix_ = pickle.load(open(f"{params['pathAnalysis']}/{params['nameMatrix']}", 'rb'))
    list_sample_ = pd.read_csv(f"{params['pathAnalysis']}/pruned_chr21.fam", sep = ' ', header = None)
    list_sample_ = list_sample_[[0]]
    info_ = []
    for row_ in range(matrix_.shape[0]):
        for col_ in range(matrix_.shape[1]):
            row_subject = list_sample_.values[row_][0]
            col_subject = list_sample_.values[col_][0]
            info_.append([row_,col_,row_subject, col_subject, matrix_[row_,col_]])
    df_ = pd.DataFrame(info_, columns = ['r_idx','c_idx','r_subject','c_subject','value'])
    df_.loc[:,'label_'] = params['label_']
    # assigning axes
    list_sample_.rename(columns = {0: 'subject_id'}, inplace = True)
    id_data_ordered = pd.merge(list_sample_, id_data, on = ['subject_id'])
    id_data_ordered.sort_values(['pop','subject_id'], inplace = True)
    id_data_ordered.reset_index(inplace = True)
    id_data_ordered.loc[:, 'desired_idx'] = id_data_ordered.index + 1
    infos_matrix = pd.merge(df_, id_data_ordered[['subject_id','pop','desired_idx']], left_on = ['r_subject'], right_on=['subject_id'])
    infos_matrix.rename(columns = {'pop':'r_pop','desired_idx':'r_desired_idx'}, inplace = True)
    infos_matrix = pd.merge(infos_matrix, id_data_ordered[['subject_id','pop','desired_idx']], left_on = ['c_subject'], right_on=['subject_id'])
    infos_matrix.rename(columns = {'pop':'c_pop','desired_idx':'c_desired_idx'}, inplace = True)
    return infos_matrix


# Setting paths
conf_paths = ["../../conf/local"]
conf_loader = ConfigLoader(conf_paths)
parameters = conf_loader.get("paths*", "paths*/**")
pathTemp = parameters['pathTemp']


paths_ = {
            'eur_orig': 'snpsParams_maf_0.05_hwe_1e-07_vif_1.5_sampParams_FIN_GBR_TSI_CEU_sex_null_lab_null_outliers_null_genes_HLA-A'
            , 'eur_univar_out': 'snpsParams_maf_0.05_hwe_1e-07_vif_1.5_sampParams_FIN_GBR_TSI_CEU_sex_null_lab_null_outliers_0.01_genes_HLA-A'
            , 'eur_multivar_out': 'snpsParams_maf_0.05_hwe_1e-07_vif_1.5_sampParams_FIN_GBR_TSI_CEU_sex_null_lab_null_outliers_0.01_genes_HLA-A_HLA-B_HLA-C_HLA-DQA1_HLA-DPB1_HLA-DRB1_HLA-DPA1_HLA-DQB1_HLA-DRA'
            , 'all_orig': 'snpsParams_maf_0.05_hwe_1e-07_vif_1.5_sampParams_FIN_GBR_TSI_CEU_YRI_sex_null_lab_null_outliers_null_genes_HLA-A'
            , 'all_univar_out': 'snpsParams_maf_0.05_hwe_1e-07_vif_1.5_sampParams_FIN_GBR_TSI_CEU_YRI_sex_null_lab_null_outliers_0.01_genes_HLA-A'
            , 'all_multivar_out': 'snpsParams_maf_0.05_hwe_1e-07_vif_1.5_sampParams_FIN_GBR_TSI_CEU_YRI_sex_null_lab_null_outliers_0.01_genes_HLA-A_HLA-B_HLA-C_HLA-DQA1_HLA-DPB1_HLA-DRB1_HLA-DPA1_HLA-DQB1_HLA-DRA'
        }
matrix_name_ = 'K_C_-1_pipelineFullChrSer.pkl'

phen_data = pd.read_csv('../../data/01_raw/hla_expression.tsv', sep = '\t')
id_data = phen_data.loc[phen_data.gene_name == 'HLA-A'].copy()
matrices_ = []
for key_ in paths_:
    folder_ = paths_[key_]
    params_it_ = {
                'pathAnalysis': f'''{pathTemp}/{folder_}'''
                , 'nameMatrix': matrix_name_
                , 'label_': key_
    }
    df_matrix = matrix_to_df(params_it_, id_data)
    matrices_.append(df_matrix)


concat_matrices = pd.concat(matrices_, axis = 0, ignore_index=True)

infos_matrix = concat_matrices.copy()
for pop_ in infos_matrix.r_pop.unique():
    print(pop_)
    print(infos_matrix.loc[(infos_matrix.r_pop == pop_)&(infos_matrix.c_pop != pop_)&(infos_matrix.r_desired_idx != infos_matrix.c_desired_idx), 'value'].max())
    print(infos_matrix.loc[(infos_matrix.r_pop == pop_)&(infos_matrix.c_pop != pop_)&(infos_matrix.r_desired_idx != infos_matrix.c_desired_idx), 'value'].min())

values_cut = infos_matrix.groupby(['label_', 'r_pop']).agg({'r_desired_idx':["min","max"]}).reset_index()
values_cut.columns = ['label','pop','min_','max_']
values_cut.loc[:, 'mean'] = (values_cut.min_ + values_cut.max_)/2

plt_matrix = pd.merge(infos_matrix, values_cut, left_on=['label_','r_pop'], right_on=['label','pop'])

theme_set(theme_bw())
p = ggplot() +\
    geom_tile(plt_matrix, aes(x='r_desired_idx',y='c_desired_idx',fill = 'value')) +\
    facet_wrap(['label_']) +\
    geom_text(plt_matrix, aes('mean', 'mean' , label = 'pop')) +\
    labs(x='Row index', y='Column index',fill = 'Value')
    
    


p.save(filename = 'matrices_by_sample.png', height=5, width=5, units = 'in', dpi=1000)

def order_comparison(x,y):
    elements_ = set([x,y])
    ordered_ = sorted(elements_)
    if len(ordered_) > 1:
        return " - ".join([ordered_[0],ordered_[1]])
    else:
        return ordered_[0]


off_diag = plt_matrix.loc[(plt_matrix.r_desired_idx != plt_matrix.c_desired_idx)].copy()
off_diag.loc[:, 'comparison_'] = off_diag.apply(lambda x: order_comparison(x['r_subject'], x['c_subject']) ,axis = 1)
off_diag.loc[:, 'comparison_pops_'] = off_diag.apply(lambda x: order_comparison(x['r_pop'], x['c_pop']) ,axis = 1)
off_diag.drop_duplicates(['comparison_', 'label_'], inplace = True)
off_diag.loc[off_diag.c_pop == off_diag.r_pop,'Case'] = "Same"
off_diag.loc[off_diag.c_pop != off_diag.r_pop,'Case'] = "Different"
p = ggplot() +\
    labs(x='Pop', y='Values',fill = 'Pop') +\
    geom_boxplot(off_diag.loc[off_diag.Case == 'Different'], aes(x = 'comparison_pops_', y = 'value',fill = 'comparison_pops_')) +\
    labs(x = "Population pairs", y = 'Value', fill = None) +\
    facet_wrap(['label_']) +\
    guides(fill=None) +\
    theme(
        axis_text_x=element_text(rotation=90, hjust=1)
)
    # geom_text(values_cut, aes(x = 'mean', y = 'mean', label = 'pop'), colour = 'white')
p.save(filename = 'comparison_values_different_by_sample.png', height=5, width=5, units = 'in', dpi=1000)

p = ggplot() +\
    labs(x='Pop', y='Values',fill = 'Pop') +\
    geom_boxplot(off_diag.loc[off_diag.Case != 'Different'], aes(x = 'comparison_pops_', y = 'value',fill = 'comparison_pops_')) +\
    labs(x = "Population pairs", y = 'Value', fill = None) +\
    facet_wrap(['label_']) +\
    guides(fill=None) +\
    theme(
        axis_text_x=element_text(rotation=90, hjust=1)
)
    # geom_text(values_cut, aes(x = 'mean', y = 'mean', label = 'pop'), colour = 'white')
p.save(filename = 'comparison_values_same_by_sample.png', height=5, width=5, units = 'in', dpi=1000)

tab_ = off_diag.groupby(['Case','comparison_pops_','label_'])['value'].describe()
tab_.to_csv('summary_comparisons_tab_by_sample.txt', sep = '|', decimal = ',')

main_diag = plt_matrix.loc[(plt_matrix.r_desired_idx == plt_matrix.c_desired_idx)].copy()
p = ggplot() +\
    labs(x='Pop', y='Values',fill = 'Pop') +\
    geom_boxplot(main_diag, aes(x = 'r_pop', y = 'value',fill = 'r_pop')) +\
    labs(x = "Population", y = 'Value', fill = None) +\
    guides(fill=None) +\
    facet_wrap(['label_']) +\
    theme(
        axis_text_x=element_text(rotation=90, hjust=1)
)
    # geom_text(values_cut, aes(x = 'mean', y = 'mean', label = 'pop'), colour = 'white')
p.save(filename = 'comparison_pops_diagonal_by_sample.png', height=5, width=5, units = 'in', dpi=1000)
