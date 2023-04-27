from kedro.config import ConfigLoader
from plotnine import *
import itertools
import numpy as np
import os
import pandas as pd
import pickle

conf_paths = ["../../conf/local"]
conf_loader = ConfigLoader(conf_paths)
parameters = conf_loader.get("paths*", "paths*/**")
temp_folder = parameters['pathTemp']


mafs = [0.05,.1]
# mafs = [0.01]
pops = [['FIN','GBR','TSI','CEU','YRI'], ['FIN','GBR','TSI','CEU']]
# pops = [['FIN','GBR','TSI','CEU','YRI']]
vifs = [1.5,2]
hwes = [1e-7,1e-4]
# hwes = [1e-7]
sexs = ['null']
labs = ['null']
genes = [
        ["HLA-A", "HLA-B", "HLA-C", "HLA-DQA1", "HLA-DPB1", "HLA-DRB1", "HLA-DPA1", "HLA-DQB1", "HLA-DRA"]
        , ['HLA-A']
        ]

outliers = [0.01, 'null']
matrices_ = {'K_C':[-1]}
sets_ = {'pipelineFullChrSer':[x for x in range(1,23)]}


dfs_ = []
file_check_ = []
for element in itertools.product(mafs , pops , vifs , hwes , sexs , labs , outliers,genes):
    maf_ = str(element[0])
    pop_ = element[1]
    vif_ = str(element[2])
    hwe_ = str(element[3])
    sex_ = element[4]
    lab_ = element[5]
    outliers_ = str(element[6])
    genes_ = element[7]
    popStr_ = "_".join(pop_)
    sexStr_ = "_".join(sex_) if sex_ != 'null' else sex_
    labStr_ = "_".join(lab_) if lab_ != 'null' else lab_
    geneStr_ = "_".join(genes_)
    folder_name_ = f'''snpsParams_maf_{maf_}_hwe_{hwe_}_vif_{vif_}_sampParams_{popStr_}_sex_{sexStr_}_lab_{labStr_}_outliers_{outliers_}_genes_{geneStr_}'''
    for file_ in os.listdir(f'''{temp_folder}/{folder_name_}'''):
        if 'resultfit_' in file_:
            df_ = pd.read_csv(f'''{temp_folder}/{folder_name_}/{file_}''', sep = ' ')
            df_.loc[:, 'fileName'] = file_
            df_.loc[:,'vif'] = vif_
            df_.loc[:,'maf'] = maf_
            df_.loc[:,'hwe'] = hwe_
            df_.loc[:,'outliers'] = outliers_
            df_.loc[:,'pop'] = 'European' if len(pop_) == 4 else 'All'
            df_.loc[:,'genes'] = 'Solo' if len(genes_) == 1 else 'All'
            df_.loc[:,'gene'] = file_.split('_')[-1].replace('.txt','')
            df_.loc[:,'type_'] = 'std' if 'std_' in file_ else 'original'
            dfs_.append(df_)
            file_check_.append(file_)


# Debug:
pd.DataFrame(file_check_, columns = ['file']).file.value_counts()


final_df = pd.concat(dfs_, axis = 0, ignore_index = True)
eigen_info_df = pickle.load(open('eigen_info.pkl','rb'))
eigen_info_df = eigen_info_df.drop(columns=['eigen']).drop_duplicates()
gene_interest = final_df.loc[final_df.gene == 'HLA-A']
df_plt = pd.merge(
        gene_interest
        , eigen_info_df
        , on = ['vif','maf','hwe','outliers','pop', 'genes','type_']
        )

df_plt.sort_values(['vif','maf','hwe','outliers'], inplace = True)
df_plt.loc[:, 'combination_'] = df_plt.apply(lambda x: f'''{x['vif']}_{x['maf']}_{x['hwe']}_{x['outliers']}''', axis = 1)
df_plt.loc[:, 'upper_bound'] = df_plt.apply(lambda x: min(x['h1'] + x['sd_h1'], 1), axis = 1)
df_plt.loc[:, 'lower_bound'] = df_plt.apply(lambda x: max(x['h1'] - x['sd_h1'], 0), axis = 1)
p_ = ggplot(df_plt.loc[df_plt.type_ == 'original'], aes(x = 'combination_', y='h1', colour = 'status', ymax = 'upper_bound', ymin = 'lower_bound')) +\
    geom_point() +\
    geom_errorbar(width=.2, position=position_dodge(.9)) +\
    facet_grid(['genes','pop'], scales='free_y') +\
    theme(figure_size=(12,9), axis_text_x=element_text(rotation=90, hjust=1)) +\
    coord_flip()

p_.save('herit_estimates')

p_ = ggplot(df_plt.loc[df_plt.type_ == 'std'], aes(x = 'combination_', y='h1', colour = 'status', ymax = 'upper_bound', ymin = 'lower_bound')) +\
    geom_point() +\
    geom_errorbar(width=.2, position=position_dodge(.9)) +\
    facet_grid(['genes','pop'], scales='free_y') +\
    theme(figure_size=(12,9), axis_text_x=element_text(rotation=90, hjust=1)) +\
    coord_flip()

p_.save('herit_estimates_std')
