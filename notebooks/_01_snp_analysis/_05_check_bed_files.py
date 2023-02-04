import pandas as pd
from multiprocessing import Pool
from kedro.config import ConfigLoader
from bed_reader import open_bed
import numpy as np

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
    totalBadColumns = len(np.where(stds_ == 0)[0])
    return totalBadColumns


def howManyBadCols(params):
    label_ = params['label_']
    info_ = []
    for chr_ in range(1,23):
        info_.append([label_, chr_, badColsComputer(params,chr_)])
    return pd.DataFrame(info_, columns = ['label', 'chr', 'badCols'])


# Setting paths
conf_paths = ["../../conf/local"]
conf_loader = ConfigLoader(conf_paths)
parameters = conf_loader.get("paths*", "paths*/**")
name_ = 'test_'
pathTemp = parameters['pathTemp']
pathTempFiles = f'''{pathTemp}/{name_}'''

originalFilesParams = {
                        'pathAnalysis': pathTempFiles,
                        'filePatt_': 'stdPruned',
                        'label_': 'Original',
                        }

dedupFilesParams = {
                        'pathAnalysis': pathTempFiles,
                        'filePatt_': 'dedupPruned',
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
tab_.to_csv('totalBadCols', index = False, sep = '|')

print(f'''
Considering a bad column is a column without variation (at all), we have for each set:

Original: 
    - average: {int(totalBadCols.loc[totalBadCols['label'] == originalFilesParams['label_'], 'badCols'].mean())} 
    - total: {int(totalBadCols.loc[totalBadCols['label'] == originalFilesParams['label_'], 'badCols'].sum())}

New set:
    - average: {int(totalBadCols.loc[totalBadCols['label'] == dedupFilesParams['label_'], 'badCols'].mean())} 
    - total: {int(totalBadCols.loc[totalBadCols['label'] == dedupFilesParams['label_'], 'badCols'].sum())}
    ''')
