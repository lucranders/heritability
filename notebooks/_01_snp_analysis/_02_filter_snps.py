import subprocess
import pandas as pd
import os
import logging
from typing import Any, Dict
import shutil
from kedro.config import ConfigLoader


# Reading dataframe
pathPhenotypeData = '../../data/01_raw/hla_expression.tsv'
pathGenotypeSample = '../../data/01_raw/SampleTot.txt'

dfPhen_ = pd.read_csv(pathPhenotypeData, sep = "\t")
samplePhen_ = dfPhen_.drop_duplicates(['subject_id'])

print(f'''
SAMPLE WITH PHENOTYPE DATA:
{samplePhen_['pop'].value_counts()}

SIZE: {samplePhen_.shape[0]}
''')

sampleGen_ = pd.read_csv(pathGenotypeSample, sep = "\t")
print(f'''
SAMPLE WITH GENOTYPE DATA:
{sampleGen_.shape[0]}
''')

sampleAnalysis = pd.merge(
                        samplePhen_,
                        sampleGen_,
                        how = 'inner'
)
print(f'''
SAMPLE TO ANALYSIS:
{sampleAnalysis['pop'].value_counts()}

TOTAL SIZE:
{sampleAnalysis.shape[0]}


LOST: 
{samplePhen_.shape[0] - sampleAnalysis.shape[0]}

{samplePhen_['pop'].value_counts() - sampleAnalysis['pop'].value_counts()}
''')


def createTempFolder(name_: str):
    # Create temporary folder
    conf_paths = ["../../conf/local"]
    conf_loader = ConfigLoader(conf_paths)
    parameters = conf_loader.get("paths*", "paths*/**")
    pathTemp = parameters['pathTemp']
    pathTempFiles = f'''{pathTemp}/{name_}'''
    proc_ = subprocess.Popen(['mkdir',pathTempFiles])
    proc_.wait()
    logging.info('Created temporary folder')
    return pathTempFiles

def selectSample(
    data: pd.DataFrame , 
    genoData: pd.DataFrame, 
    sampParams: dict, 
    pathTempFiles: str
    ) -> Dict[str, Any]:
    """
    Node for selecting the desired sample to work with.
    The parameters are taken from conf/project/parameters.yml.
    The data and the parameters will be loaded and provided to this function
    automatically when the pipeline is executed and it is time to run this node.
    """
    if sampParams['pop'] != None:
        data = data.loc[data['pop'].isin(sampParams['pop'])]
    if sampParams['sex'] != None:
        data = data.loc[data['sex'].isin(sampParams['sex'])]
    if sampParams['lab'] != None:
        data = data.loc[data['lab'].isin(sampParams['lab'])]
    # Merged data - intersection between genotyped and phenotyped individuals
    dfMerge = genoData.merge(data, how = 'inner')
    sizeOriginalDf = dfMerge.drop_duplicates(['subject_id']).shape[0]
    # logging.info('The original sample size is: ' + str(sizeOriginalDf))
    # Saving sample on temp folder
    if sampParams['outliers'] != None:
        from scipy.linalg import inv
        import numpy as np
        from scipy.stats import chi2
        cut_ = 1 - sampParams['outliers']
        gene_names = dfMerge.gene_name.unique()
        horizontalDf = dfMerge.pivot(index=['subject_id','sex','lab','pop'], columns=['gene_name'], values='tpm').reset_index()
        centeredValues = horizontalDf[gene_names] - horizontalDf[gene_names].mean()
        centeredValues.index = horizontalDf.subject_id
        covMatrix_ = np.cov(horizontalDf[gene_names].values.T)
        invCovMatrix_ = inv(covMatrix_)
        invCovMatrix_ = pd.DataFrame(invCovMatrix_)
        invCovMatrix_.index = gene_names
        invCovMatrix_.columns = gene_names
        mahalanobisMatrix = centeredValues @ invCovMatrix_ @ centeredValues.T
        mahalanobisFinal = np.diag(mahalanobisMatrix)
        mahalanobisFinal = pd.DataFrame(mahalanobisFinal)
        mahalanobisFinal.loc[:,'subject_id'] = horizontalDf.subject_id
        mahalanobisFinal.rename(columns = {0:'MahalanobisD'},inplace = True)
        mahalanobisFinal.loc[:,'Percentile'] = mahalanobisFinal.MahalanobisD.apply(lambda x: chi2.cdf(x,(len(gene_names))))
        mahalanobisFinal.loc[:,'cut'] = mahalanobisFinal.Percentile.apply(lambda x: 1 if x >= cut_ else 0)
        dfSample = mahalanobisFinal.loc[mahalanobisFinal.cut == 0 ,['subject_id']].copy().drop_duplicates()
    else:
        dfSample = dfMerge.loc[:,['subject_id']].copy().drop_duplicates()
        mahalanobisFinal = dfMerge
    print(f'The new sample size is: {dfSample.shape[0]}')
    dfSample.loc[:,'1'] = dfSample.subject_id
    dfSample.to_csv(pathTempFiles + '/sample.txt', index = False, header= False, sep = ' ')
    final_ = {'selectedSample':dfSample.loc[:,['subject_id']] , 'sizeOriginalDf':sizeOriginalDf , 'sizeFinalDf':dfSample.shape[0] , 'originalDf': dfMerge, 'MahalanobisDf': mahalanobisFinal,'pathAnalysis': pathTempFiles}
    return final_

sampParams = {'pop': None,'sex':None,'lab':None,'outliers':None}
pathTemp = createTempFolder('test_')

sample = selectSample(
                        dfPhen_.loc[dfPhen_.subject_id.isin(sampleAnalysis.subject_id)],
                        sampleGen_,
                        sampParams,
                        pathTemp
                        )
print(f'''
SIZE OF ORIGINAL SAMPLE:
{sample['sizeOriginalDf']}
SIZE OF NEW SAMPLE
{sample['sizeFinalDf']}
''')

def checkLogSizes(path_,listNames_):
    listSizes = []
    for name_ in listNames_:
        try:
            size = os.path.getsize(path_ + '/' + name_)
        except:
            size = 0
    
        listSizes.append(size > 0)
    return listSizes



def filterSnps(pathTempFiles:str, snpsParams: dict, nameFiles:str) -> None:
    conf_paths = ["../../conf/local"]
    conf_loader = ConfigLoader(conf_paths)
    parameters = conf_loader.get("paths*", "paths*/**")
    pathPlink = parameters['plink']
    pathVcf = parameters['pathVcfFiles']
    query_1 = pathPlink + " --vcf $input"
    if snpsParams['vif'] != None:
        query_1 += " --indep 50 5 " + str(snpsParams['vif'])
    if snpsParams['maf'] != None:
        query_1 += " --maf " + str(snpsParams['maf'])
    if snpsParams['hwe'] != None:
        query_1 += " --hwe " + str(snpsParams['hwe'])
    query_1 += f" --keep {pathTempFiles}/sample.txt --mind 0.05 --geno 0.05 --vcf-half-call missing --out $filtered"
    # check if files were already created
    checkCreatedFiles_ = checkLogSizes(pathTempFiles,[f'list_filt_snps_chr{chr_}_{nameFiles}.prune.in' for chr_ in range(1,23)])
    cond = sum(checkCreatedFiles_) < 22
    if cond:
        # Create .bed files - 22 chromossomes
        for chr_ in range(1,23):
            sbatchFile_1 = f'''input={pathVcf}/ALL.chr{chr_}_GRCh38.genotypes.20170504.vcf.gz
    filtered={pathTempFiles}/list_filt_snps_chr{chr_}_{nameFiles}
    eval {query_1}
            '''
            with open(pathTempFiles + '/filterSnps' + str(chr_) + '.sh', 'w') as f:
                f.write(sbatchFile_1)
                f.write('\n')
            cmd = ['sh', pathTempFiles + '/filterSnps' + str(chr_) + '.sh']
            print(sbatchFile_1)
            p = subprocess.Popen(cmd)
            p.wait()
    return None


vif_ = 1.5
maf_ = .05
hwe_ = 10**-4

filterSnps(pathTemp,{'vif': vif_, 'maf': maf_, 'hwe': hwe_},'test_')

for file_ in os.listdir(pathTemp):
    if '.prune.in' in file_:
        shutil.copyfile(f'{pathTemp}/{file_}', file_)
