import pandas as pd
from kedro.config import ConfigLoader
from bed_reader import open_bed
import numpy as np
from sklearn.preprocessing import StandardScaler
from plotnine import *
from sklearn.impute import SimpleImputer

def howManyVariants(params):
    pathAnalysis = params['pathAnalysis']
    filePatt_ = params['filePatt_']
    label_ = params['label_']
    info_ = []
    for chr_ in range(1,23):
        f_ = open(f'{pathAnalysis}/{filePatt_}_chr{chr_}.log', 'r')
        for line_ in f_:
            # Pattern to look for how many variants were selected
            if 'variants and' in line_:
                variants_ = int(line_.split(' ')[0])
                info_.append([label_, chr_, variants_])
    return pd.DataFrame(info_, columns = ['label', 'chr', 'qt'])

def read_bed_make_sample(params, df_count, chr_, alpha_):
    bedFile_ = f'''{params['pathAnalysis']}/{params['filePatt_']}_chr{chr_}.bed'''
    bed = open_bed(bedFile_)
    # Reading matrix
    val = bed.read()
    imp_mean = SimpleImputer(missing_values=np.nan, strategy='mean')
    imp_mean.fit(val)
    tempMatrix_ = imp_mean.transform(val)
    ones_ = np.ones([tempMatrix_.shape[0],1]) @ np.ones([1,tempMatrix_.shape[0]])
    centered_matrix = tempMatrix_ - ((1/tempMatrix_.shape[0]) * np.matmul(ones_,tempMatrix_))
    orig_stds_ = tempMatrix_.std(axis = 0)
    stds_ = orig_stds_ ** alpha_
    D_alpha = np.diag(stds_)
    scaled_matrix = centered_matrix @ D_alpha
    val_scaled = pd.DataFrame(scaled_matrix)
    val_scaled.loc[:,'idx'] = val_scaled.index
    melt_data = pd.melt(val_scaled, id_vars = ['idx'])
    chr_count_ = df_count.loc[df_count.chr == chr_]
    n_ = chr_count_.samp.values
    sum_idx = chr_count_.sum_idx.values[0]
    sample_vars = [x for x in melt_data.variable.drop_duplicates().sample(n_, random_state = 20230318)]
    melt_data = melt_data.loc[melt_data.variable.isin(sample_vars)]
    dict_ = {}
    count_ = 0
    for i_ in range(len(sample_vars)):
        dict_[sample_vars[i_]] = count_
        count_ += 1
    melt_data.loc[:, 'idx_var_plt'] = melt_data.variable.apply(lambda x: dict_[x] + sum_idx)
    return melt_data, orig_stds_, stds_

# Setting paths
conf_paths = ["../../conf/local"]
conf_loader = ConfigLoader(conf_paths)
parameters = conf_loader.get("paths*", "paths*/**")
pathTemp = parameters['pathTemp']

params = {
                        'pathAnalysis': f'''{pathTemp}/test_0_1''',
                        'filePatt_': 'pruned',
                        'label_': 'Without duplicates',
                        }

frac_ = .01

count_selected_variants = howManyVariants(params)
count_selected_variants.loc[:, 'rel'] = count_selected_variants.qt / count_selected_variants.qt.sum()
tot_ = count_selected_variants.qt.sum()
n_ = int(frac_ * tot_)
count_selected_variants.loc[:, 'samp'] = count_selected_variants.rel * n_
count_selected_variants.loc[:, 'samp'] = count_selected_variants.loc[:, 'samp'].astype(int)
count_selected_variants.loc[:, 'shift_idx'] = count_selected_variants.samp.shift(1)
count_selected_variants.loc[:, 'sum_idx'] = count_selected_variants.shift_idx.cumsum()
count_selected_variants.fillna(0, inplace = True)

sample_ = pd.read_csv(f'''{params["pathAnalysis"]}/sample.txt''', sep = ' ', header = None, usecols = [0])
sample_.columns = ['subject_id']
sample_.loc[:, 'desired_idx'] = sample_.index + 1
pop_info = pd.read_csv('../../data/01_raw/hla_expression.tsv', sep = '\t')
merge_ = pd.merge(pop_info[['subject_id', 'pop']].drop_duplicates(), sample_)
merge_.loc[:, 'idx'] = merge_.index

dict_stds_ = {}
for alpha_ in [-2,-1,0,1,2]:
    dfs_ = []
    dict_stds_temp_ = {}
    for chr_ in range(1,23):
        read_ = read_bed_make_sample(params, count_selected_variants, chr_, alpha_)
        dfs_.append(read_[0])
        dict_stds_temp_[chr_] = {'std_orig': read_[1], 'std_transformed': read_[2]}
    df_concat_ = pd.concat(dfs_, axis = 0, ignore_index = True)
    dict_stds_[alpha_] = dict_stds_temp_
    df_concat_merge_ = pd.merge(df_concat_, merge_)
    df_concat_merge_.sort_values(['desired_idx'], inplace = True)
    values_cut = df_concat_merge_.groupby(['pop']).agg({'desired_idx':['min','max']}).reset_index()
    values_cut.columns = ['pop','min_','max_']
    values_cut.loc[:, 'mean'] = (values_cut.min_ + values_cut.max_)/2
    p = ggplot() +\
        geom_tile(df_concat_merge_, aes(x = 'idx_var_plt', y = 'desired_idx',fill = 'value')) +\
        labs(x = "SNPs", y = 'Index', fill = "Value") +\
        geom_text(values_cut, aes(x = 0, y = 'mean', label = 'pop'), colour = 'white', alpha = .75) +\
        geom_vline(count_selected_variants, aes(xintercept = 'sum_idx'), colour = 'red', size = .05, alpha = .75) +\
        geom_hline(values_cut, aes(yintercept = 0), colour = 'red', size = .05, alpha = .75) +\
        geom_hline(values_cut, aes(yintercept = 'max_'), colour = 'red', size = .05, alpha = .75) +\
        theme(
            axis_text_x = element_blank(),

    )
    p.save(filename = f'matrix_alpha_{alpha_}_2.png', height=5, width=5, units = 'in', dpi=1000)