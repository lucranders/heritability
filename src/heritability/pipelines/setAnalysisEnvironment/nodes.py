from typing import Any, Dict
import subprocess
import os
import time
from datetime import datetime
import pandas as pd
import numpy as np
import pandas as pd
import logging
from kedro.config import ConfigLoader

logging.info('Current working directory:' + os. getcwd())
def createTempFolder(snpsParams: dict):
    # Create temporary folder
    conf_paths = ["conf/local"]
    conf_loader = ConfigLoader(conf_paths)
    parameters = conf_loader.get("paths*", "paths*/**")
    pathTemp = parameters['pathTemp'] 
    pathTempFiles = pathTemp + '/TempBed_maf_' + str(snpsParams['maf']) + "_hwe_" + str(snpsParams['hwe']) + "_vif_" + str(snpsParams['vif'])
    proc_ = subprocess.Popen(['mkdir',pathTempFiles])
    proc_.wait()
    logging.info('Created temporary folder')
    return pathTempFiles

def selectSample(data: pd.DataFrame , genoData: pd.DataFrame, sampParams: dict, pathTempFiles:str) -> Dict[str, Any]:
    """Node for selecting the desired sample to work with.
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
    logging.info('The original sample size is: ' + str(sizeOriginalDf))
    # Saving sample on temp folder
    if sampParams['outliers'] != None:
        from scipy.linalg import inv
        import numpy as np
        from scipy.stats import chi2
        cut_ = 1 - sampParams['outliers']
        gene_names = data.gene_name.unique()
        horizontalDf = data.pivot(index=['subject_id','sex','lab','pop'], columns=['gene_name'], values='tpm').reset_index()
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
    sizeFinalDf = dfSample.shape[0]
    logging.info('The new sample size is: ' + str(sizeFinalDf))
    dfSample.loc[:,'1'] = dfSample.subject_id
    dfSample.to_csv(pathTempFiles + '/sample.txt', index = False, header= False, sep = ' ')
    final_ = {'selectedSample':dfSample.loc[:,['subject_id']] , 'sizeOriginalDf':sizeOriginalDf , 'sizeFinalDf':sizeFinalDf , 'originalDf': dfMerge}
    return final_
def checkLogSizes(path_,name_,ext_):
    listSizes = []
    for chr_ in range(1,23):
        try:
            size = os.path.getsize(path_ + '/' + name_ + str(chr_) + ext_)
        except:
            size = 0
    
        listSizes.append(size > 0)
    return listSizes
# Update status about process
def updateLog(path_,status_,file_):
    if status_ == 0:
        msg_ = str(datetime.now()) + ": Waiting..."
    else:
        msg_ = str(datetime.now()) + ": Done!"
    f = open(path_ + '/' + file_, "a")
    f.write(msg_ + '\n')
    f.close()
def monitoringProcess(pathTempFiles,name_,ext_,statusFile_):
    # Check if all 22 .bed files are created
    sizes = checkLogSizes(pathTempFiles,name_,ext_)
    cond = sum(sizes) < 22
    # If not, waits until so
    while cond:
        sizes = checkLogSizes(pathTempFiles,name_,ext_)
        cond = sum(sizes) < 22
        if cond:
            # Updates status (inside tmp folder)
            updateLog(pathTempFiles,0,statusFile_ + ".txt")
            time.sleep(10)
        else:
            # Updates status (inside tmp folder)
            updateLog(pathTempFiles,1,statusFile_ + ".txt")
def createBedFiles(pathTempFiles: str , snpsParams: dict, selectedSample: dict) -> None:
    conf_paths = ["conf/local"]
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
    # query_ += " --keep " + pathTempFiles + '/sample.txt' + ' --mind 0.05 --geno 0.05 --vcf-half-call missing --make-bed --out 
    query_1 += " --keep " + pathTempFiles + '/sample.txt' + ' --mind 0.05 --geno 0.05 --vcf-half-call missing --out $filtered --noweb' 
    # Create .bed files - 22 chromossomes
    sbatchFile_1 = f'''#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --partition=short
#SBATCH --array=1-22
#SBATCH --time=480
#SBATCH --mem=16gb
#SBATCH --job-name=filterSnps
#SBATCH --nodelist=darwin
#SBATCH --output=/scratch/genevol/users/lucas/out/filterSnps.txt

input={pathVcf}/ALL.chr"$SLURM_ARRAY_TASK_ID"_GRCh38.genotypes.20170504.vcf.gz
filtered={pathTempFiles}/list_filt_snps_chr"$SLURM_ARRAY_TASK_ID"
eval {query_1}
    '''
    with open(pathTempFiles + '/filterSnps.sbatch', 'w') as f:
        f.write(sbatchFile_1)
        f.write('\n')
    cmd = ['sbatch', pathTempFiles + '/filterSnps.sbatch']
    # subprocess.Popen(cmd)
    # monitoringProcess(pathTempFiles,'list_filt_snps_chr','.prune.in','filteredSnpsStatus')
    query_2 = pathPlink + " --vcf $input --vcf-half-call missing --extract $filtered --make-bed --out $bed --noweb"
    sbatchFile_2 = f'''#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --partition=short
#SBATCH --array=1-22
#SBATCH --time=480
#SBATCH --mem=16gb
#SBATCH --job-name=bedMaker
#SBATCH --nodelist=darwin
#SBATCH --output=/scratch/genevol/users/lucas/out/bedMaker.txt

input={pathVcf}/ALL.chr"$SLURM_ARRAY_TASK_ID"_GRCh38.genotypes.20170504.vcf.gz
filtered={pathTempFiles}/list_filt_snps_chr"$SLURM_ARRAY_TASK_ID".prune.in
bed={pathTempFiles}/chr"$SLURM_ARRAY_TASK_ID"
eval {query_2}
    '''
    with open(pathTempFiles + '/bedMaker.sbatch', 'w') as f:
        f.write(sbatchFile_2)
        f.write('\n')
    cmd = ['sbatch', pathTempFiles + '/bedMaker.sbatch']
    # subprocess.Popen(cmd)   
    monitoringProcess(pathTempFiles,'chr','.bed','bedStatus')
    return 1