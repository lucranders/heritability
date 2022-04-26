from typing import Any, Dict
import subprocess
import os
import time
from datetime import datetime
import pandas as pd
import numpy as np
import pandas as pd

def createTempFolder(snpsParams: dict, pathTemp: str):
    # Create temporary folder
    pathTempFiles = pathTemp + '/TempBed_maf_' + str(snpsParams['maf']) + "_hwe_" + str(snpsParams['hwe']) + "_vif_" + str(snpsParams['vif'])
    proc_ = subprocess.Popen(['mkdir',pathTempFiles])
    proc_.wait()
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
        mahalanobisFinal.loc[:,'Percentile'] = mahalanobisFinal.MahalanobisD.apply(lambda x: chi2.cdf(x,(len(gene_names)-1)))
        mahalanobisFinal.loc[:,'cut'] = mahalanobisFinal.Percentile.apply(lambda x: 1 if x >= cut_ else 0)
        dfSample = mahalanobisFinal.loc[mahalanobisFinal.cut == 0 ,['subject_id']].drop_duplicates()
    else:
        dfSample = dfMerge.loc[:,['subject_id']].drop_duplicates()
    dfSample.loc[:,'1'] = dfSample.subject_id
    dfSample.to_csv(pathTempFiles + '/sample.txt', index = False, header= False, sep = ' ')
    return dfMerge
def checkLogSizes(path_):
    listSizes = []
    for chr_ in range(1,23):
        try:
            size = os.path.getsize(path_ + '/chr' + str(chr_) + '.log')
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
def monitoringSnpSelectionBedFiles(pathTempFiles):
    # Check if all 22 .bed files are created
    sizes = checkLogSizes(pathTempFiles)
    cond = sum(sizes) < 22
    # If not, waits until so
    while cond:
        sizes = checkLogSizes(pathTempFiles)
        cond = sum(sizes) < 22
        if cond:
            # Updates status (inside tmp folder)
            updateLog(pathTempFiles,0,"bedStatus.txt")
            time.sleep(10)
        else:
            # Updates status (inside tmp folder)
            updateLog(pathTempFiles,1,"bedStatus.txt")
def createBedFiles(pathPlink:str , pathTempFiles: str , pathVcf: str , snpsParams: dict, selectedSample: pd.DataFrame) -> None:
    query_ = pathPlink + " --vcf $input"
    if snpsParams['vif'] != None:
        query_ += " --indep 50 5 " + str(snpsParams['vif'])
    if snpsParams['maf'] != None:
        query_ += " --maf " + str(snpsParams['maf'])
    if snpsParams['hwe'] != None:
        query_ += " --hwe " + str(snpsParams['hwe'])
    # query_ += " --keep " + pathTempFiles + '/sample.txt' + ' --mind 0.05 --geno 0.05 --vcf-half-call missing --make-bed --out 
    query_ += " --keep " + pathTempFiles + '/sample.txt' + ' --mind 0.05 --geno 0.05 --vcf-half-call missing --make-bed --out $bed' 
    # Create .bed files - 22 chromossomes
    cmd = ['qsub', '-v' ,'query_=' + query_ + ',tempPath=' + pathTempFiles + ',pathVcf=' + pathVcf , 'filterSnps.sh']
    subprocess.Popen(cmd)
    monitoringSnpSelectionBedFiles(pathTempFiles)
    return 1