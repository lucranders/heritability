import pandas as pd
from kedro.config import ConfigLoader
from bed_reader import open_bed
import numpy as np
import pickle

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


def badColsComputer(params, chr_):
    pathAnalysis = params['pathAnalysis']
    filePatt_ = params['filePatt_']
    bedFile_ = f'{pathAnalysis}/{filePatt_}_chr{chr_}.bed'
    bed = open_bed(bedFile_)
    # Reading matrix
    val = bed.read()
    # Computing std, ignoring nan values
    stds_ = np.nanstd(val, axis = 0)
    noVar = np.where(stds_ == 0)[0]
    totalBadColumns = len(noVar)
    return totalBadColumns, noVar


def howManyBadCols(params):
    label_ = params['label_']
    info_ = []
    for chr_ in range(1,23):
        info_.append([label_, chr_, badColsComputer(params,chr_)[0]])
    return pd.DataFrame(info_, columns = ['label', 'chr', 'badCols'])

def whichBadCols(params):
    label_ = params['label_']
    info_ = []
    for chr_ in range(1,23):
        info_.append([label_, chr_, badColsComputer(params,chr_)[1]])
    return pd.DataFrame(info_, columns = ['label', 'chr', 'badCols'])


# Setting paths
conf_paths = ["../../conf/local"]
conf_loader = ConfigLoader(conf_paths)
parameters = conf_loader.get("paths*", "paths*/**")
pathTemp = parameters['pathTemp']

originalFilesParams = {
                        'pathAnalysis': f'''{pathTemp}/test_0_0''',
                        'filePatt_': 'pruned',
                        'label_': 'Original',
                        }

dedupFilesParams = {
                        'pathAnalysis': f'''{pathTemp}/test_0_1''',
                        'filePatt_': 'pruned',
                        'label_': 'Without duplicates',
                        }

selectedVariants = pd.concat([howManyVariants(originalFilesParams), howManyVariants(dedupFilesParams)], ignore_index= True, axis = 0)
tab_ = pd.pivot(selectedVariants, index = 'chr',columns='label', values = 'qt').reset_index()
tab_.loc[:, 'diff'] = tab_[originalFilesParams['label_']] - tab_[dedupFilesParams['label_']]
print(f'''
On average, {int(tab_['diff'].mean())} variants are dropped
Considering all chromosomes, the total is {int(tab_['diff'].sum())}
''')

totalBadCols = pd.concat([howManyBadCols(originalFilesParams), howManyBadCols(dedupFilesParams)], ignore_index= True, axis = 0)
tab_ = pd.pivot(totalBadCols, index = 'chr',columns='label', values = 'badCols').reset_index()
tab_.to_csv(f'''{originalFilesParams['pathAnalysis']}/totalBadCols.txt''', index = False, sep = '|')

print(f'''
Considering a bad column is a column without variation (at all), we have for each set:

Original: 
    - average: {int(totalBadCols.loc[totalBadCols['label'] == originalFilesParams['label_'], 'badCols'].mean())} 
    - total: {int(totalBadCols.loc[totalBadCols['label'] == originalFilesParams['label_'], 'badCols'].sum())}

New set:
    - average: {int(totalBadCols.loc[totalBadCols['label'] == dedupFilesParams['label_'], 'badCols'].mean())} 
    - total: {int(totalBadCols.loc[totalBadCols['label'] == dedupFilesParams['label_'], 'badCols'].sum())}
    ''')

index_bad_cols_orig = whichBadCols(originalFilesParams)
pickle.dump(index_bad_cols_orig, open(f'''{originalFilesParams['pathAnalysis']}/index_bad_cols.pkl''','wb'))
index_bad_cols_orig.to_csv(f'''{originalFilesParams['pathAnalysis']}/index_bad_cols.txt''', index = False, sep = '|')