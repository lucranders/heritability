import pandas as pd
import os
import shutil

mergedFiles_ = []
filescwd = [x for x in os.listdir()]
for chr_ in range(1,23):
    prunedSnpsFile = [x for x in filescwd if '.prune.in' in x and f'chr{chr_}_' in x][0]
    snpsInfoFile = [x for x in filescwd if f't{chr_}.txt' in x][0]
    print(f'''
    Pruned snps file: {prunedSnpsFile}
    Snps info file: {snpsInfoFile}
    ''')
    # Read pruned files and rename columns
    auxPruned = pd.read_csv(prunedSnpsFile, sep = '\t', header = None)
    auxPruned.rename(columns = {0: 'SNP'}, inplace = True)
    auxPruned.loc[:, 'CHR'] = chr_
    auxPruned.loc[:, 'aux_'] = 1
    # Count how many duplicates the snp has in the file
    auxPruned2 = auxPruned.groupby(['SNP','CHR']).aux_.sum().reset_index()
    auxPruned2.rename(columns = {'aux_': 'duplicatesPruned'}, inplace = True)
    print('''
    TOTAL OF DUPLICATED SNPS ALREADY COMPUTED
    ''')
    # Read snps info and rename columns
    auxSnpsInfo = pd.read_csv(snpsInfoFile, sep = '\t', header = None)
    auxSnpsInfo.rename(columns = {0: 'CHR', 1: 'POS', 2: 'SNP'}, inplace = True)
    auxSnpsInfo = auxSnpsInfo.loc[auxSnpsInfo.SNP.isin(auxPruned2.SNP)]
    print('''
    SNPS OF INTEREST ALREADY FILTERED (INFO)
    ''')
    auxSnpsInfo.loc[:, 'aux_'] = 1
    # Count how many duplicates the snp has in the file (by snp)
    auxSnpsInfo2 = auxSnpsInfo.groupby(['SNP','CHR']).aux_.sum().reset_index()
    auxSnpsInfo2.rename(columns = {'aux_': 'duplicatesInfoBySNP'}, inplace = True)
    print('''
    TOTAL OF DUPLICATED SNPS INFO BY SNP ALREADY COMPUTED
    ''')
    # Count how many duplicates the snp has in the file (by snp and pos)
    auxSnpsInfo3 = auxSnpsInfo.groupby(['POS','SNP','CHR']).aux_.sum().reset_index()
    auxSnpsInfo3.rename(columns = {'aux_': 'duplicatesInfoBySNPPos'}, inplace = True)
    print('''
    TOTAL OF DUPLICATED SNPS INFO BY SNP AND POSITION ALREADY COMPUTED
    ''')
    # Merge infos
    auxSnpsInfo4 = pd.merge(
                            auxSnpsInfo3
                            , auxSnpsInfo2
                            , how = 'left'
                            )
    # merge infos and pruned snps
    merge_ = pd.merge(
                    auxPruned2
                    , auxSnpsInfo4
                    , how = 'left'
                    )
    mergedFiles_.append(merge_)


compiledInfo = pd.concat(mergedFiles_, axis = 0, ignore_index = True)

tab_ = compiledInfo.groupby(['CHR','duplicatesInfoBySNPPos', 'duplicatesInfoBySNP']).SNP.count().reset_index()
print(tab_.loc[tab_.duplicatesInfoBySNPPos == tab_.duplicatesInfoBySNP].shape[0] == tab_.shape[0])
# Conclusion: The same snp is not located in two different places - check
del(tab_)
tab_ = compiledInfo.groupby(['CHR','duplicatesInfoBySNPPos', 'duplicatesPruned']).SNP.count().reset_index()
print(tab_.loc[tab_.duplicatesPruned > tab_.duplicatesInfoBySNPPos].shape[0] == 0)
# Conclusion: When there are duplicates, it does not happen to be more duplicates on selection than on available information - check
print(round(100 * tab_.loc[tab_.duplicatesInfoBySNPPos == 1].SNP.sum() / compiledInfo.shape[0], 2))
# If we drop all the duplicate snps, we will keep 99.86% of all snps


conf_paths = ["../../conf/local"]
conf_loader = ConfigLoader(conf_paths)
parameters = conf_loader.get("paths*", "paths*/**")
pathTemp = parameters['pathTemp']
for chr_ in range(1,23):
    # aux_ = compiledInfo.loc[(compiledInfo.duplicatesInfoBySNPPos == 1)&(compiledInfo.CHR == chr_)]
    file_ = f'new_list_filt_snps_{chr_}.prune.in'
    # aux_.loc[:,['SNP']].to_csv(file_, index = False, header = None, sep = '\t')
    shutil.copyfile(file_ , f'{pathTemp}/test_/{file_}')
