import pickle
import itertools
import os
import numpy as np
import pandas as pd
from kedro.config import ConfigLoader
from plotnine import *

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
        if '.pkl' in file_:
            if 'correction_non_negative' not in file_:
                print(genes_, len(genes_))
                matrix_ = pickle.load(open(f'{temp_folder}/{folder_name_}/{file_}','rb'))
                w,v = np.linalg.eig(matrix_)
                df_ = pd.DataFrame(w,columns=['eigen'])
                df_.loc[:, 'fileName'] = file_
                df_.loc[:,'vif'] = vif_
                df_.loc[:,'maf'] = maf_
                df_.loc[:,'hwe'] = hwe_
                df_.loc[:,'outliers'] = outliers_
                df_.loc[:,'pop'] = 'European' if len(pop_) == 4 else 'All'
                df_.loc[:,'genes'] = 'Solo' if len(genes_) == 1 else 'All'
                df_.loc[:, 'n_sample'] = df_.shape[0]
                if df_.eigen.min() < 0:
                    df_.loc[:,'status'] = 'negative'
                else:
                    df_.loc[:,'status'] = 'positive'
                dfs_.append(df_)


finalDf_ = pd.concat(dfs_, axis = 0, ignore_index = True)
finalDf_.sort_values(['vif','maf','hwe','outliers'], inplace = True)
finalDf_.loc[:, 'combination_'] = finalDf_.apply(lambda x: f'''{x['vif']}_{x['maf']}_{x['hwe']}_{x['outliers']}''', axis = 1)

p_ = ggplot(finalDf_, aes(x = 'combination_', y='eigen', colour = 'status')) +\
    geom_boxplot() +\
    facet_wrap(['genes','pop'], scales='free_y') +\
    theme(figure_size=(12,9), axis_text_x=element_text(rotation=90, hjust=1))

p_.save('boxplot_eigen_values')

