from typing import Any, Dict
import subprocess
import os
import time
from datetime import datetime
from pathlib import Path
import pandas as pd
import numpy as np
import pandas as pd


def createTmpFolder(pathTemp: str, snpsParams: dict):
    # Create temp folder
    path_ = pathTemp + '/TempBed_maf_' + str(snpsParams['maf']) + "_hwe_" + str(snpsParams['hwe']) + "_vif_" + str(snpsParams['vif'])
    subprocess.Popen(['mkdir',path_])
    return path_

def selectSample(data: pd.DataFrame, sampParams: dict, pathTempFiles: str) -> Dict[str, Any]:
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
    # Merge with genotypes samples
    # 
    # 
    # 
    # Temporary
    # Merged data - intersection between genotyped and phenotyped individuals
    dfMerge = data.copy()
    # Saving sample on temp folder
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
def createBedFiles(pathPlink:str , pathTempFiles: str , pathVcf: str , snpsParams: dict):
    print(pathPlink,pathTempFiles,pathVcf,snpsParams)
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
    # cmd = ['qsub', '-v' ,'query_=' + query_ + ',tempPath=' + pathTempFiles , 'filterSnps.sh']
    cmd = ['sh', '-v' ,'query_=' + query_ + ',tempPath=' + pathTempFiles + ',pathVcf=' + pathVcf + ',chr=' + str(22) , 'filterSnps.sh']
    subprocess.Popen(cmd)

def monitoringSnpSelection(pathTempFiles):
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

# Calculate ZZ' for the given set of snps
def calculateGCTA(self,nameFile = None,nameMatrix = None):
    if nameFile == None:
        refChrs = 'chrs'
    if nameMatrix != None:
        self.nameMatrix_ = 'GCTA_' + nameMatrix
    else:
        self.nameMatrix_ = 'GCTA'
    if self.sample_ != None:
        cmd = [self.pathGCTA_ + '/gcta64', '--mbfile' ,self.path_ + '/' + refChrs + '.txt','--keep',self.path_ + '/sample.txt','--make-grm','--out',self.path_+'/' + self.nameMatrix_,'--thread-num',self.threads_]
    else:
        cmd = [self.pathGCTA_ + '/gcta64', '--mbfile' ,self.path_ + '/' + refChrs + '.txt','--make-grm','--out',self.path_+'/' + self.nameMatrix_,'--thread-num',self.threads_]
    subprocess.Popen(cmd)
    # check whether GRM binaries already exists - means process is finished
    check_ = Path(self.path_ + '/' + self.nameMatrix_ + '.grm.bin')
    cond = check_.is_file()
    while not cond:
        cond = check_.is_file()
        if not cond:
            # Updates status (inside tmp folder)
            updateLog(self.path_,0,'GRM' + self.nameMatrix_ + 'Status.txt')
            time.sleep(10)
        else:
            # Updates status (inside tmp folder)
            updateLog(self.path_,1,'GRM' + self.nameMatrix_ + 'Status.txt')

