import pickle
import itertools
import os
import numpy as np
import pandas as pd

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
        ]

outliers = [0.01, 'null']
formulas = [{'fixed': '~pop + lab', 'random': "null"}]
matrices_ = {'K_C':[-1]}
sets_ = {'pipelineFullChrSer':[x for x in range(1,23)]}


configs_ = []
for element in itertools.product(mafs , pops , vifs , hwes , sexs , labs , outliers,genes,formulas):
    maf_ = str(element[0])
    pop_ = element[1]
    vif_ = str(element[2])
    hwe_ = str(element[3])
    sex_ = element[4]
    lab_ = element[5]
    outliers_ = str(element[6])
    genes_ = element[7]
    formula_ = element[8]
    popStr_ = "_".join(pop_)
    sexStr_ = "_".join(sex_) if sex_ != 'null' else sex_
    labStr_ = "_".join(lab_) if lab_ != 'null' else lab_
    geneStr_ = "_".join(genes_)
    folder_name_ = f'''snpsParams_maf_{maf_}_hwe_{hwe_}_vif_{vif_}_sampParams_{popStr_}_sex_{sexStr_}_lab_{labStr_}_outliers_{outliers_}_genes_{geneStr_}'''
    configs_.append([
        maf_
        , vif_
        , hwe_
        , outliers_
        , popStr_
        , folder_name_
    ])

configs_df_ = pd.DataFrame(configs_, columns = ['maf','vif','hwe','outliers','pop', 'folder'])

def generateData(root_):
    dfs_ = []
    for file_ in os.listdir(root_):
        if '.pkl' in file_:
            matrix_ = pickle.load(open(f'{root_}/{file_}','rb'))
            nameFile = file_.split('raw_')[1]
            qtSnps = int(nameFile.replace('.pkl',''))
            nameFile = file_.split('config_')[1]
            configNumber_ = nameFile[0]
            w,v = np.linalg.eig(matrix_)
            df_ = pd.DataFrame(w,columns=['eigen'])
            df_.loc[:, 'qtSnps'] = qtSnps
            config_ = configs_[f'config_{configNumber_}']
            df_.loc[:, 'configNumber'] = configNumber_
            df_.loc[:, 'fileName'] = file_
            df_.loc[:,'vif'] = config_['vif']
            df_.loc[:,'maf'] = config_['maf']
            df_.loc[:,'hwe'] = config_['hwe']
            df_.loc[:,'config'] = 'hwe:' + df_.hwe.astype(str) + ' maf:' + df_.maf.astype(str) + ' vif:' + df_.vif.astype(str) + ' snps:' + df_.qtSnps.astype(str)
            if df_.eigen.min() < 0:
                df_.loc[:,'negative'] = 1
            else:
                df_.loc[:,'negative'] = 0
            dfs_.append(df_)
    finalDf_ = pd.concat(dfs_, axis = 0, ignore_index = True)
    return finalDf_    

def bxPltKC(root_, ax_, title_):
    dfBxPlt = generateData(root_)
    sns.boxplot(data = dfBxPlt, y = 'eigen', x = 'config', ax = ax_, )
    ax_.tick_params(labelrotation=90)
    ax_.set_title(title_)
    # medians = dfBxPlt.groupby(['config'])['eigen'].median().values
    numberSnps = dfBxPlt.groupby(['config'])['qtSnps'].max().values.astype(str)
    negativeFlag = dfBxPlt.groupby(['config'])['negative'].max().values
    mins_ = dfBxPlt.groupby(['config'])['eigen'].min().values
    numberSnps = ["snps: " + i for i in numberSnps]
    pos = range(len(numberSnps))
    for tick,label in zip(pos,ax_.get_xticklabels()):
        # ax_.text(pos[tick],
        #         medians[tick] + .03,
        #         numberSnps[tick],
        #         horizontalalignment='center',
        #         size='x-small',
        #         color='red',
        #         rotation = 'vertical',
        #         weight='semibold')
        if negativeFlag[tick] == 1:
            ax_.text(pos[tick],
                    mins_[tick] + 0.03,
                    np.round(mins_[tick],7),
                    horizontalalignment='center',
                    size='x-small',
                    color='red',
                    weight='semibold')
    return ax_
    
fig, (ax0,ax1,ax2,ax3) = plt.subplots(1, 4,  sharex = False,figsize=(16,16))
bxPltKC('/scratch/genevol/users/lucas/sample1/', ax0, 'Everyone')
bxPltKC('/scratch/genevol/users/lucas/sample2/', ax1, 'Everyone but outliers')
bxPltKC('/scratch/genevol/users/lucas/sample3/', ax2, 'Europeans')
bxPltKC('/scratch/genevol/users/lucas/sample4/', ax3, 'Europeans but outliers')
# * 1 - Everyone, no outlier filters;
# * 2 - Everyone, filtering out outliers based on Mahalanobis distance (percentile 99%);
# * 3 - Europeans, no outlier filters;
# * 4 - Europeans, filtering out outliers based on Mahalanobis distance (percentile 99%)
fig.tight_layout()