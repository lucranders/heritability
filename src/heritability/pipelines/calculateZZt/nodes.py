"""
This is a boilerplate pipeline 'calculateZZt'
generated using Kedro 0.17.6
"""
import subprocess
from pathlib import Path
import time
from datetime import datetime
import numpy as np
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
import os
from kedro.config import ConfigLoader
from os.path import exists

def updateLog(path_,status_,file_):
    if status_ == 0:
        msg_ = str(datetime.now()) + ": Waiting..."
    else:
        msg_ = str(datetime.now()) + ": Done!"
    f = open(path_ + '/' + file_, "a")
    f.write(msg_ + '\n')
    f.close()
def createChrRef(pathAnalysis,listChrs = None,nameFile = 'chrs'):
    # Create references to calculate GCTA
    if listChrs == None:
        for chr_ in range(1,23):
            f = open(pathAnalysis + '/' + nameFile +'.txt', "a")
            f.write(pathAnalysis+'/chr' + str(chr_) + '\n')
            f.close()
    else:
        for chr_ in listChrs:
            f = open(pathAnalysis + '/' + nameFile +'.txt', "a")
            f.write(pathAnalysis+'/chr' + str(chr_) + '\n')
            f.close()

# Calculate ZZ' for the given set of snps
def calculateGCTA(matrixes:dict,selectedSample: dict ,numThreads:int, check_bed_files_exists:int):
    pathAnalysis = selectedSample['pathAnalysis']
    print(pathAnalysis)
    for nameMatrixIt_ in matrixes:
        nameMatrix = 'GCTA_' + nameMatrixIt_
        if not exists(pathAnalysis+'/' + nameMatrix + '.grm.bin'):
            conf_paths = ["conf/local"]
            conf_loader = ConfigLoader(conf_paths)
            parameters = conf_loader.get("paths*", "paths*/**")
            pathGCTA = parameters['gcta']
            listChrs = matrixes[nameMatrixIt_]
            nameFile = 'chromosomes' + nameMatrixIt_
            # print(pathAnalysis,listChrs,nameFile)
            createChrRef(pathAnalysis=pathAnalysis,listChrs=listChrs,nameFile=nameFile)
            # print('ok')
            cmd = [pathGCTA, '--mbfile' ,pathAnalysis + '/' + nameFile + '.txt','--keep',pathAnalysis + '/sample.txt','--make-grm','--out',pathAnalysis+'/' + nameMatrix,'--thread-num',str(numThreads)]
            # print(cmd)
            subprocess.Popen(cmd)
            # check whether GRM binaries already exists - means process is finished
            # print(pathAnalysis + '/' + nameMatrix + '.grm.bin')
            check_ = Path(pathAnalysis + '/' + nameMatrix + '.grm.bin')
            # print(check_)
            cond = check_.is_file()
            while not cond:
                cond = check_.is_file()
                if not cond:
                    # Updates status (inside tmp folder)
                    updateLog(pathAnalysis,0,'GRM' + nameMatrix + 'Status.txt')
                    time.sleep(10)
                else:
                    # Updates status (inside tmp folder)
                    updateLog(pathAnalysis,1,'GRM' + nameMatrix + 'Status.txt')
        else:
            print("ZZt already calculated!")
    return 1
def correctGRM(selectedSample: dict, matrixes: dict, check_GCTA_calculated:int):
    pathAnalysis = selectedSample['pathAnalysis']
    dictGenes = dict()
    for nameMatrixIt_ in matrixes:
        nameMatrix = 'GCTA_' + nameMatrixIt_
        # Reads and stores as a variable the calculated matrix
        pandas2ri.activate()
        readRDS = robjects.r['readRDS']
        if not exists(pathAnalysis + '/' + nameMatrix + '_correction.rds'):
            call_ = subprocess.Popen(['Rscript','src/heritability/pipelines/calculateZZt/matrixCorrection.r',pathAnalysis,nameMatrix])
            call_.wait()
            correctedGCTA = readRDS(pathAnalysis + '/' + nameMatrix + '_correction.rds')
        else:
            print('ZZt already corrected!')
            correctedGCTA = readRDS(pathAnalysis + '/' + nameMatrix + '_correction.rds')
        dictGenes[nameMatrixIt_] = correctedGCTA
    return dictGenes
