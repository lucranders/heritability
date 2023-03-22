from kedro.config import ConfigLoader
import numpy as np
import pandas as pd
import struct
import subprocess
from plotnine import *

def calculateGCTA(args_):
    f'''
    Builds the GCTA matrix.
    args_: dict
        pathRef_: path to files (reference to bed files + samples)
        nameFile: name of the reference file (bed files)
        nameMatrix: name of the final matrix (output)
        numThreads_: number of threads to calculate matrix

    '''
    pathAnalysis = args_['pathAnalysis']
    filePatt_ = args_['filePatt_']
    nameMatrix_ = args_['nameMatrix']
    numThreads_ = 10
    conf_paths = ["../../conf/local"]
    conf_loader = ConfigLoader(conf_paths)
    parameters = conf_loader.get("paths*", "paths*/**")
    pathGCTA_ = parameters['gcta']
    # Write ref file
    file_ = open(f'{pathAnalysis}/REF_BEDS_{filePatt_}', 'w')
    [file_.write(f'{pathAnalysis}/{filePatt_}_chr{chr_}\n') for chr_ in range(1,23)]
    file_.close()
    # command to build matrix
    cmd = [pathGCTA_, '--mbfile' ,f'{pathAnalysis}/REF_BEDS_{filePatt_}', '--keep', f'{pathAnalysis}/sample.txt' ,'--make-grm','--out',f'{pathAnalysis}/{nameMatrix_}','--thread-num',str(numThreads_)]
    # Execute command and wait for the end
    p = subprocess.Popen(cmd)
    p.wait()
    return "Done!"

# Read GRM as a matrix
def readBinaryGCTA(file_):
    t_ = pd.read_csv(f'{file_}.grm.id', sep = '\t', header = None)
    n_ = t_.shape[0]
    biteSize = 4
    totalElements = int(n_ * (n_ + 1) / 2)
    values_ = []
    with open(f'{file_}.grm.bin','rb') as binFile_:
        for i_ in range(totalElements):
            values_.append(struct.unpack('f',binFile_.read(biteSize))[0])
    dfIndex = pd.DataFrame([x for x in range(n_)], columns = ['index_'])
    dfIndex.loc[:,'index_'] = dfIndex.index_ + 1
    dfIndex.loc[:, 'indexC_'] = dfIndex.index_.cumsum()
    dfIndex.loc[:, 'indexDiagonal'] = dfIndex.indexC_
    indexDiagonal = dfIndex.indexDiagonal.values - 1
    diagonalMatrix = [values_[x] for x in indexDiagonal]
    offDiagonalMatrix = [values_[x] for x in range(len(values_)) if x not in indexDiagonal]
    # Build matrix from binary files
    gctaMatrix_ = np.zeros([n_,n_])
    cont_ = 0
    contDiag_ = 0
    for col_ in range(n_):
        for row_ in range(n_):
            if row_ < col_:
                gctaMatrix_[row_,col_] = offDiagonalMatrix[cont_]
                gctaMatrix_[col_,row_] = offDiagonalMatrix[cont_]
                cont_ += 1
            if row_ == col_:
                gctaMatrix_[row_,col_] = diagonalMatrix[contDiag_]
                contDiag_ += 1
    return gctaMatrix_

# Extract eigen values
def eigenValues(params):
    matrix_ = readBinaryGCTA(f"{params['pathAnalysis']}/{params['nameMatrix']}")
    label_ = params['label_']
    w,v = np.linalg.eig(matrix_)
    del v
    df_ = pd.DataFrame(w,columns=[label_])
    return df_

# Setting paths
conf_paths = ["../../conf/local"]
conf_loader = ConfigLoader(conf_paths)
parameters = conf_loader.get("paths*", "paths*/**")
name_ = 'test_0'
pathTemp = parameters['pathTemp']
pathTempFiles = f'''{pathTemp}/{name_}'''

originalFilesParams = {
                        'pathAnalysis': pathTempFiles,
                        'filePatt_': 'stdPruned',
                        'nameMatrix': 'original',
                        'label_': 'Original'
                        }

dedupFilesParams = {
                        'pathAnalysis': pathTempFiles,
                        'filePatt_': 'dedupPruned',
                        'nameMatrix': 'withoutDuplicates',
                        'label_': 'Without duplicates'
                        }
# Calculate matrix
calculateGCTA(originalFilesParams)
calculateGCTA(dedupFilesParams)

# Compile eigen values
ordered_eigen_original_set = eigenValues(originalFilesParams).sort_values(originalFilesParams['label_'], ascending = False)
ordered_eigen_original_set.loc[:, 'idx'] = range(1,ordered_eigen_original_set.shape[0] + 1)
ordered_eigen_dedup_set = eigenValues(dedupFilesParams).sort_values(dedupFilesParams['label_'], ascending = False)
ordered_eigen_dedup_set.loc[:, 'idx'] = range(1,ordered_eigen_dedup_set.shape[0] + 1)
dfEigen = pd.merge(ordered_eigen_original_set, ordered_eigen_dedup_set, how = 'inner', on = ['idx'])
dfEigen.loc[:,'dif'] = dfEigen[originalFilesParams['label_']] - dfEigen[dedupFilesParams['label_']]
dfEigen.loc[:, 'negative1'] = dfEigen[originalFilesParams['label_']].apply(lambda x: '<0' if x < 0 else '>= 0')
dfEigen.loc[:, 'negative2'] = dfEigen[dedupFilesParams['label_']].apply(lambda x: '<0' if x < 0 else '>= 0')

p = ggplot(dfEigen, aes(x=originalFilesParams['label_'],y=dedupFilesParams['label_'],colour = 'negative1', shape = 'negative2')) +\
    labs(x='x', y='y') +\
    geom_point(size=0.1) +\
    geom_abline(intercept = 0, slope = 1, alpha = .5)+\
    labs(colour = 'Original set', shape = 'New set', x = 'Eigen values (original set)', y = 'Eigen values (new set)')
    
p.save(filename = 'eigenValues.png', height=5, width=5, units = 'in', dpi=1000)

theme_set(theme_bw())
p = ggplot(dfEigen, aes(x='idx',y='dif',colour = 'negative1', shape = 'negative2')) +\
    labs(x='x', y='y') +\
    geom_point(size=0.1) +\
    labs(colour = 'Original set', shape = 'New set', x = 'Index of eigenvalue', y = 'Difference between eigenvalues')
    

p.save(filename = 'eigenValuesDif.png', height=5, width=5, units = 'in', dpi=1000)


def matrix_to_df(params):
    matrix_ = readBinaryGCTA(f"{params['pathAnalysis']}/{params['nameMatrix']}")
    list_sample_ = pd.read_csv(f"{params['pathAnalysis']}/{params['nameMatrix']}.grm.id", sep = '\t', header = None)
    list_sample_.drop(columns = [1], inplace = True)
    list_sample_ = list_sample_.values
    info_ = []
    for row_ in range(matrix_.shape[0]):
        for col_ in range(matrix_.shape[1]):
            # row_subject = list_sample_.iloc[row_]['subject_id']
            # col_subject = list_sample_.iloc[col_]['subject_id']
            row_subject = list_sample_[row_][0]
            col_subject = list_sample_[col_][0]
            info_.append([row_,col_,row_subject, col_subject, matrix_[row_,col_]])
    df_ = pd.DataFrame(info_, columns = ['r_idx','c_idx','r_subject','c_subject','value'])
    df_.loc[:,'label_'] = params['label_']
    return df_

matrix_values_orig = matrix_to_df(originalFilesParams)
matrix_values_dedup = matrix_to_df(dedupFilesParams)
concat_matrices = pd.concat([matrix_values_orig, matrix_values_dedup], axis = 0, ignore_index=True)

phen_data = pd.read_csv('../../data/01_raw/hla_expression.tsv', sep = '\t')
id_data = phen_data.loc[phen_data.gene_name == 'HLA-A'].copy()
sample_data = pd.read_csv(f'{pathTempFiles}/sample.txt', sep =' ', header = None)
sample_data.drop(columns = [1], inplace = True)
sample_data.rename(columns = {0: 'subject_id'}, inplace = True)
id_data_ordered = pd.merge(sample_data, id_data, on = ['subject_id'])
id_data_ordered.loc[:, 'desired_idx'] = id_data_ordered.index + 1
# id_data_ordered.loc[:, 'aux'] = 1
# id_data_ordered.loc[:, 'order_by_pop'] = id_data_ordered.groupby(['pop']).aux.cumsum()
# max_groups_id_ = id_data_ordered.groupby(['pop']).order_by_pop.max().reset_index()
# max_groups_id_.loc[:,'freq'] = max_groups_id_.order_by_pop.cumsum()
# max_groups_id_.loc[:, 'bounds'] = max_groups_id_.freq.shift(1)
# max_groups_id_.fillna(0, inplace = True)
# id_data_ordered.loc[:, 'desired_idx'] = id_data_ordered.order_by_pop.apply(lambda x: x + max_groups_id_)
infos_matrix = pd.merge(concat_matrices, id_data_ordered[['subject_id','pop','desired_idx']], left_on = ['r_subject'], right_on=['subject_id'])
infos_matrix.rename(columns = {'pop':'r_pop','desired_idx':'r_desired_idx'}, inplace = True)
infos_matrix = pd.merge(infos_matrix, id_data_ordered[['subject_id','pop','desired_idx']], left_on = ['c_subject'], right_on=['subject_id'])
infos_matrix.rename(columns = {'pop':'c_pop','desired_idx':'c_desired_idx'}, inplace = True)

infos_matrix.sort_values([])



for pop_ in infos_matrix.r_pop.unique():
    print(pop_)
    print(infos_matrix.loc[(infos_matrix.r_pop == pop_)&(infos_matrix.c_pop != pop_)&(infos_matrix.r_desired_idx != infos_matrix.c_desired_idx), 'value'].max())
    print(infos_matrix.loc[(infos_matrix.r_pop == pop_)&(infos_matrix.c_pop != pop_)&(infos_matrix.r_desired_idx != infos_matrix.c_desired_idx), 'value'].min())

values_cut = infos_matrix.groupby(['r_pop']).agg({'r_desired_idx':["max","min"]}).reset_index()
values_cut.columns = ['pop','max_','min_']
values_cut.loc[:, 'mean'] = (values_cut.min_ + values_cut.max_)/2

theme_set(theme_bw())
p = ggplot() +\
    labs(x='Row index', y='Column index',fill = 'Value') +\
    geom_tile(infos_matrix, aes(x='r_desired_idx',y='c_desired_idx',fill = 'value')) +\
    facet_wrap(['label_']) +\
    geom_text(values_cut, aes(x = 'mean', y = 'mean', label = 'pop'), colour = 'white')


p.save(filename = 'matrices_.png', height=5, width=5, units = 'in', dpi=1000)

def order_comparison(x,y):
    elements_ = set([x,y])
    ordered_ = sorted(elements_)
    if len(ordered_) > 1:
        return " - ".join([ordered_[0],ordered_[1]])
    else:
        return ordered_[0]


off_diag = infos_matrix.loc[(infos_matrix.r_desired_idx != infos_matrix.c_desired_idx)&(infos_matrix.label_ == 'Without duplicates')].copy()
off_diag.loc[:, 'comparison_'] = off_diag.apply(lambda x: order_comparison(x['r_subject'], x['c_subject']) ,axis = 1)
off_diag.loc[:, 'comparison_pops_'] = off_diag.apply(lambda x: order_comparison(x['r_pop'], x['c_pop']) ,axis = 1)
off_diag.drop_duplicates(['comparison_'], inplace = True)
off_diag.loc[off_diag.c_pop == off_diag.r_pop,'Case'] = "Same"
off_diag.loc[off_diag.c_pop != off_diag.r_pop,'Case'] = "Different"
p = ggplot() +\
    labs(x='Pop', y='Values',fill = 'Pop') +\
    geom_boxplot(off_diag.loc[off_diag.Case == 'Different'], aes(x = 'comparison_pops_', y = 'value',fill = 'comparison_pops_')) +\
    labs(x = "Population pairs", y = 'Value', fill = None) +\
    guides(fill=None) +\
    theme(
        axis_text_x=element_text(rotation=90, hjust=1)
)
    # geom_text(values_cut, aes(x = 'mean', y = 'mean', label = 'pop'), colour = 'white')
p.save(filename = 'comparison_values_different.png', height=5, width=5, units = 'in', dpi=1000)

p = ggplot() +\
    labs(x='Pop', y='Values',fill = 'Pop') +\
    geom_boxplot(off_diag.loc[off_diag.Case != 'Different'], aes(x = 'comparison_pops_', y = 'value',fill = 'comparison_pops_')) +\
    labs(x = "Population pairs", y = 'Value', fill = None) +\
    guides(fill=None) +\
    theme(
        axis_text_x=element_text(rotation=90, hjust=1)
)
    # geom_text(values_cut, aes(x = 'mean', y = 'mean', label = 'pop'), colour = 'white')
p.save(filename = 'comparison_values_same.png', height=5, width=5, units = 'in', dpi=1000)

tab_ = off_diag.groupby(['Case','comparison_pops_'])['value'].describe()
tab_.to_csv('summary_comparisons_tab', sep = '|', decimal = ',')

main_diag = infos_matrix.loc[(infos_matrix.r_desired_idx == infos_matrix.c_desired_idx)&(infos_matrix.label_ == 'Without duplicates')].copy()
p = ggplot() +\
    labs(x='Pop', y='Values',fill = 'Pop') +\
    geom_boxplot(main_diag, aes(x = 'r_pop', y = 'value',fill = 'r_pop')) +\
    labs(x = "Population", y = 'Value', fill = None) +\
    guides(fill=None) +\
    theme(
        axis_text_x=element_text(rotation=90, hjust=1)
)
    # geom_text(values_cut, aes(x = 'mean', y = 'mean', label = 'pop'), colour = 'white')
p.save(filename = 'comparison_pops_diagonal.png', height=5, width=5, units = 'in', dpi=1000)
