"""
This is a boilerplate pipeline 'calculateZZt'
generated using Kedro 0.17.6
"""
import time
import subprocess
import rpy2.robjects as robjects
import pickle
import pandas as pd
import os
import numpy as np
from pathlib import Path
from os.path import exists
from multiprocessing import Pool
from kedro.config import ConfigLoader
from datetime import datetime
from bed_reader import open_bed

class calculateMatrixes:
    def __init__(self,params):
        self.chromosomes_ = params["chromosomes"]
        self.pathAnalysis_ = params["pathAnalysis"] + '/'
        self.alpha = params["alpha"]
        self.cores_ = params["cores"]
        self.sample_ = pd.read_csv(self.pathAnalysis_ + 'sample.txt', header = None, sep = ' ', usecols = [0], names = ['subject_id'])
    def Bool_(self,x,comp_):
        if x in comp_:
            return True
        else:
            return False
    def readBedsMP(self,chr_):
        bed = open_bed(self.pathAnalysis_ + "chr" + str(chr_) + ".bed")
        val = bed.read()
        boolArray_ = [self.Bool_(x,self.sample_.subject_id.values) for x in bed.iid]
        desiredSample = bed.iid[boolArray_]
        filterRows = val[boolArray_]
        return chr_, filterRows
    def readBeds(self):
        matrixes_ = {}
        p = Pool(self.cores_)
        with p:
            results = p.map(self.readBedsMP,self.chromosomes_)
        for i_ in range(len(self.chromosomes_)):
            result_ = results[i_]
            chr_ = result_[0]
            mat_ = result_[1]
            matrixes_[chr_] = mat_
        self.matrixes_ = matrixes_
    def eliminateBadColumns(self):
        correctedMatrixes = {}
        for chr_ in self.chromosomes_:
            tempMatrix = self.matrixes_[chr_].copy()
            # Calculates std after subtracting values by the mean and replacing by 0 nan values 
            # (equivalent to replace by the mean)
            stds_ = np.nan_to_num(tempMatrix - tempMatrix.mean(axis = 0)).std(axis = 0)
            checkVar = np.where(stds_ == 0)
            if len(checkVar) > 0:
                correctedMatrix_ = np.nan_to_num(np.delete(tempMatrix,checkVar, axis = 1))
            else:
                correctedMatrix_ = np.nan_to_num(tempMatrix.copy())
            correctedMatrixes[chr_] = correctedMatrix_.copy()
        self.correctedMatrixDef = correctedMatrixes
    def calculateMatrix(self):
        ncolsDict = {}
        zzts_ = {}
        for i_,chr_ in enumerate(self.chromosomes_):
            tempMatrix_ = self.correctedMatrixDef[chr_].copy()
            means_ = tempMatrix_.mean(axis = 0)
            stds_ = tempMatrix_.std(axis = 0)
            ncols_ = tempMatrix_.shape[1]
            ncolsDict[chr_] = ncols_
            Z_ = ((tempMatrix_ - means_) * ((stds_)**self.alpha))
            ZZt_ = Z_ @ Z_.T
            zzts_[chr_] = ZZt_
            if i_ == 0:
                zztMatrixes = ZZt_.copy()
            else: 
                zztMatrixes += ZZt_
        self.zztMatrixes = zztMatrixes.copy()
        self.ncolsMatrixes = ncolsDict
        self.zzts_ = zzts_
    def calculateFinalMatrix(self):
        for i_,chr_ in enumerate(self.chromosomes_):
            if i_ == 0:
                nTot = self.ncolsMatrixes[chr_]
            else:
                nTot += self.ncolsMatrixes[chr_]
        self.nTot = nTot
        self.zztFinal = self.zztMatrixes / self.nTot


def calculate_K_C_alpha(matrixes:dict,selectedSample: dict ,parametersMatrix:dict, check_bed_files_exists:int):
    pathAnalysis = selectedSample['pathAnalysis']
    print(pathAnalysis)
    for alpha_ in parametersMatrix["alphas"]:
        for nameMatrixIt_ in matrixes:
            nameMatrix = 'K_C_' + str(alpha_) + '_' + nameMatrixIt_
            if not exists(pathAnalysis+'/' + nameMatrix + '.pkl'):
                params = {}
                params["pathAnalysis"] = pathAnalysis
                params["alpha"] = alpha_
                params["cores"] = parametersMatrix["cores"]
                params["chromosomes"] = matrixes[nameMatrixIt_]
                calculator = calculateMatrixes(params)
                calculator.readBeds()
                calculator.eliminateBadColumns()
                calculator.calculateMatrix()
                calculator.calculateFinalMatrix()
                desiredMatrix = calculator.zztFinal
                with open(pathAnalysis+'/' + nameMatrix + '.pkl','wb') as wf_:
                    pickle.dump(desiredMatrix, wf_)
            else:
                print('K_C_' + str(alpha_) + '_' + nameMatrixIt_ + " already calculated!")
    return 1

def correct_K_C_alpha(matrixes: dict, selectedSample: dict, parametersMatrix: dict, check_GCTA_calculated:int):
    pathAnalysis = selectedSample['pathAnalysis']
    for alpha_ in parametersMatrix["alphas"]:
        for nameMatrixIt_ in matrixes:
            nameMatrix = 'K_C_' + str(alpha_) + '_' + nameMatrixIt_
            # Reads and stores as a variable the calculated matrix
            if not exists(pathAnalysis + '/' + nameMatrix + '_correction.pkl'):
                with open(pathAnalysis+'/' + nameMatrix + '.pkl','rb') as rf_:
                    desiredMatrix = pickle.load(rf_)
                try:
                    choleskyDecomp_ = np.linalg.cholesky(desiredMatrix)
                    a = 1
                except:
                    a = 0
                if a == 0:
                    newMatrix = desiredMatrix.copy()
                    while a == 0:
                        eigenValues, eigenVectors = np.linalg.eig(newMatrix)
                        minEigen = np.min(abs(eigenValues))
                        reg_ = min(abs(minEigen), 1e-6)
                        newDiag = np.diag(eigenValues + reg_)
                        newMatrix = eigenVectors @ newDiag @ eigenVectors.T
                        try: 
                            np.linalg.cholesky(newMatrix)
                            a = 1
                        except:
                            a = 0
                    with open(pathAnalysis + '/' + nameMatrix + '_correction.pkl','wb') as wf_:
                        pickle.dump(newMatrix, wf_)
                else:
                    with open(pathAnalysis + '/' + nameMatrix + '_correction.pkl','wb') as wf_:
                        pickle.dump(desiredMatrix, wf_)
            else:
                print(nameMatrix + ' already corrected!')
    return 1