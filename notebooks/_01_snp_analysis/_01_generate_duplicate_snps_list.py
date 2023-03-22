import pandas as pd

for chr_ in range(1,23):
    aux_ = pd.read_csv(f'test_0/test{chr_}.txt', sep = '\t', header = None)
    aux_.columns = ['chr','pos','snp']
    tab_ = aux_.groupby(['pos','snp']).chr.count().reset_index()
    tab_.loc[tab_.chr > 1,['snp']].to_csv(f'test_0/duplicates_chr_{chr_}.txt', header = None, sep='\t', index = False)
    tab_.loc[tab_.chr > 1,['snp','chr']].to_csv(f'test_0/duplicates_counts__chr_{chr_}.txt', sep='\t', index = False)

def read_dups_(x):
    return pd.read_csv(f'test_0/duplicates_chr_{x}.txt', sep = '\t', header = None)

list_snps_ = pd.concat(map(read_dups_, [x for x in range(1,23)]), axis = 0, ignore_index= True)
list_snps_.to_csv('test_0/remove_snps.txt', sep='\t', index = False, header = None)
