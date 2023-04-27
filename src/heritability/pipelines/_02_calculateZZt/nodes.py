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
import struct
import numpy as np
from pathlib import Path
from os.path import exists
import itertools
from multiprocessing import Pool
from kedro.config import ConfigLoader
from datetime import datetime
from bed_reader import open_bed
from sklearn.impute import SimpleImputer
import logging

logging.basicConfig(filename='logs/journals/calculate_zzt.log', encoding='utf-8', level=logging.DEBUG)

class calculateMatrixes:
    def __init__(self,params):
        self.chromosomes_ = params["chromosomes"]
        self.pathAnalysis_ = params["pathAnalysis"] + '/'
        self.alpha = params["alpha"]
        self.sample_ = pd.read_csv(self.pathAnalysis_ + 'sample.txt', header = None, sep = ' ', usecols = [0], names = ['subject_id'])
    def Bool_(self,x,comp_):
        if x in comp_:
            return True
        else:
            return False
    def readBedsMP(self,chr_):
        bed = open_bed(self.pathAnalysis_ + "pruned_chr" + str(chr_) + ".bed")
        val = bed.read()
        logging.warning(f'''
        On chromosome {chr_}
            There are {np.where(pd.isna(val), 0, 1).sum()} filled values
            Imputing {np.where(pd.isna(val), 1, 0).sum()} missing values with the mean of each snp
        ''')
        # Imputing missing values with means, for each snp value
        imp_mean = SimpleImputer(missing_values=np.nan, strategy='mean')
        imp_mean.fit(val)
        tempMatrix_ = imp_mean.transform(val)
        ones_ = np.ones([tempMatrix_.shape[0],1]) @ np.ones([1,tempMatrix_.shape[0]])
        centered_matrix = tempMatrix_ - ((1/tempMatrix_.shape[0]) * np.matmul(ones_,tempMatrix_))
        stds_ = tempMatrix_.std(axis = 0) ** self.alpha
        D_alpha = np.diag(stds_)
        scaled_matrix = centered_matrix @ D_alpha
        ZZt_ = scaled_matrix @ scaled_matrix.transpose()
        n_snps_ = scaled_matrix.shape[1]
        logging.info(f'''
        Matrix on chromosome {chr_} is ready!
        There were {n_snps_} snps
        ''')
        # boolArray_ = [self.Bool_(x,self.sample_.subject_id.values) for x in bed.iid]
        # desiredSample = bed.iid[boolArray_]
        # filterRows = val[boolArray_]
        return chr_, ZZt_, n_snps_
    def readBeds(self):
        matrices_ = {}
        p = Pool(len(self.chromosomes_))
        logging.info(f"""
        Reading the following chromosomes and preparing to build a single matrix:
        {self.chromosomes_}
        """)
        with p:
            results = p.map(self.readBedsMP,self.chromosomes_)
        for i_ in range(len(self.chromosomes_)):
            result_ = results[i_]
            chr_ = result_[0]
            mat_ = result_[1]
            n_snps_ = result_[2]
            matrices_[chr_] = {'matrix': mat_, 'n_snps': n_snps_}
        self.matrices_ = matrices_
        logging.info(f"""
        Done! Moving to the next step.
        """)
    # def eliminateBadColumns(self):
    #     correctedMatrixes = {}
    #     for chr_ in self.chromosomes_:
    #         tempMatrix = self.matrices_[chr_].copy()
    #         # Calculates std after subtracting values by the mean and replacing by 0 nan values 
    #         # (equivalent to replace by the mean)
    #         stds_ = np.nan_to_num(tempMatrix - tempMatrix.mean(axis = 0)).std(axis = 0)
    #         checkVar = np.where(stds_ == 0)
    #         if len(checkVar) > 0:
    #             correctedMatrix_ = np.nan_to_num(np.delete(tempMatrix,checkVar, axis = 1))
    #         else:
    #             correctedMatrix_ = np.nan_to_num(tempMatrix.copy())
    #         correctedMatrixes[chr_] = correctedMatrix_.copy()
    #     self.correctedMatrixDef = correctedMatrixes
    def calculateMatrix(self):
        ncolsDict = {}
        logging.info(f'''
        Summing all matrices of  set {self.chromosomes_}
        ''')
        for i_,chr_ in enumerate(self.chromosomes_):
            if i_ == 0:
                zztMatrixes = self.matrices_[chr_]['matrix'].copy()
                nTot = self.matrices_[chr_]['n_snps']
            else: 
                zztMatrixes += self.matrices_[chr_]['matrix']
                nTot += self.matrices_[chr_]['n_snps']
            ncolsDict[chr_] = self.matrices_[chr_]['n_snps']
            logging.info(f'''
            Matrix on chromosome {chr_} is summed!
            ''')
        self.zztMatrixes = zztMatrixes.copy()
        self.ncolsMatrixes = ncolsDict
        self.nTot = nTot
        self.zztFinal = self.zztMatrixes / self.nTot
        logging.info(f'''
        The final matrix is ready!
        It is built with {self.nTot} snps
        ''')
    # def calculateFinalMatrixDiagCorrection(self):
    #     ncolsDict = {}
    #     zzts_ = {}
    #     for i_,chr_ in enumerate(self.chromosomes_):
    #         tempMatrix_ = self.correctedMatrixDef[chr_].copy()
    #         means_ = tempMatrix_.mean(axis = 0)
    #         stds_ = tempMatrix_.std(axis = 0)
    #         ncols_ = tempMatrix_.shape[1]
    #         ncolsDict[chr_] = ncols_
    #         # Calculating full ZZt
    #         Z_ = ((tempMatrix_ - means_) * ((stds_)**(-1)))
    #         ZZt_ = Z_ @ Z_.T
    #         # According to Yang 2010, there is a bias in the main diagonal
    #         # Calculating correction for matrix diagonal, which is:
    #         # Aijj = 1 + [x_{ij}^2 - (1 + 2p_{i})x_{ij} + 2p_{i}^2]/(2p_{i}(1-p_{i}))
    #         # But this result holds only for alpha = -1 (for now)
    #         tmpMatrix_ = ((tempMatrix_ - 1) * ((stds_)**(-1)))
    #         tmpMatrix2_ = Z_ @ tmpMatrix_.T
    #         # We are interested only on main diag
    #         mainDiag_ = np.diag(tmpMatrix2_)
    #         # Calculating matrix from additive effects
    #         ZZtCorrect = ZZt_.copy()
    #         # Fill the main diagonal with 0
    #         np.fill_diagonal(ZZtCorrect, 0)
    #         # Transforming diagonal as diagonal matrix
    #         mainDiagMatrix_ = np.diag(mainDiag_)
    #         # Replacing the main diagonal of matrix calculated from additive effects by correction
    #         ZZtCorrect = ZZtCorrect + mainDiagMatrix_
    #         zzts_[chr_] = ZZtCorrect
    #         if i_ == 0:
    #             zztMatrixes = ZZtCorrect.copy()
    #         else: 
    #             zztMatrixes += ZZtCorrect
    #     self.zztMatrixes = zztMatrixes.copy()
    #     self.ncolsMatrixes = ncolsDict
    #     self.zzts_ = zzts_



# def calculate_K_C_alpha(matrixes:dict,selectedSample: dict ,parametersMatrix:dict, check_bed_files_exists:int):
#     pathAnalysis = selectedSample['pathAnalysis']
#     print(pathAnalysis)
#     for alpha_ in parametersMatrix["alphas"]:
#         for nameMatrixIt_ in matrixes:
#             nameMatrix = 'K_C_' + str(alpha_) + '_' + nameMatrixIt_
#             if not exists(pathAnalysis+'/' + nameMatrix + '.pkl'):
#                 params = {}
#                 params["pathAnalysis"] = pathAnalysis
#                 params["alpha"] = alpha_
#                 params["cores"] = parametersMatrix["cores"]
#                 params["chromosomes"] = matrixes[nameMatrixIt_]
#                 calculator = calculateMatrixes(params)
#                 calculator.readBeds()
#                 calculator.eliminateBadColumns()
#                 calculator.calculateMatrix()
#                 calculator.calculateFinalMatrix()
#                 desiredMatrix = calculator.zztFinal
#                 with open(pathAnalysis+'/' + nameMatrix + '.pkl','wb') as wf_:
#                     pickle.dump(desiredMatrix, wf_)
#             else:
#                 print('K_C_' + str(alpha_) + '_' + nameMatrixIt_ + " already calculated!")
#     return 1

# Function equivalent to calculate_K_C_alpha; Only difference is the calculateMatrix step
# And matrix name
# def calculate_GCTA(matrixes:dict,selectedSample: dict ,parametersMatrix:dict, check_bed_files_exists:int):
#     pathAnalysis = selectedSample['pathAnalysis']
#     print(pathAnalysis)
#     for nameMatrixIt_ in matrixes:
#         nameMatrix = 'GCTA_' + nameMatrixIt_
#         if not exists(pathAnalysis+'/' + nameMatrix + '.pkl'):
#             params = {}
#             params["pathAnalysis"] = pathAnalysis
#             params["cores"] = parametersMatrix["cores"]
#             params["alpha"] = 'GCTA'
#             params["chromosomes"] = matrixes[nameMatrixIt_]
#             calculator = calculateMatrixes(params)
#             calculator.readBeds()
#             calculator.eliminateBadColumns()
#             calculator.calculateFinalMatrixDiagCorrection()
#             calculator.calculateFinalMatrix()
#             desiredMatrix = calculator.zztFinal
#             with open(pathAnalysis+'/' + nameMatrix + '.pkl','wb') as wf_:
#                 pickle.dump(desiredMatrix, wf_)
#         else:
#             print('GCTA_' + nameMatrixIt_ + " already calculated!")
#     return 1


# def correct_K_C_alpha(matrixes: dict, selectedSample: dict, parametersMatrix: dict, check_GCTA_calculated:int):
#     pathAnalysis = selectedSample['pathAnalysis']
#     for alpha_ in parametersMatrix["alphas"]:
#         for nameMatrixIt_ in matrixes:
#             nameMatrix = 'K_C_' + str(alpha_) + '_' + nameMatrixIt_
#             # Reads and stores as a variable the calculated matrix
#             if not exists(pathAnalysis + '/' + nameMatrix + '_correction.pkl'):
#                 with open(pathAnalysis+'/' + nameMatrix + '.pkl','rb') as rf_:
#                     desiredMatrix = pickle.load(rf_)
#                 try:
#                     choleskyDecomp_ = np.linalg.cholesky(desiredMatrix)
#                     a = 1
#                 except:
#                     a = 0
#                 if a == 0:
#                     newMatrix = desiredMatrix.copy()
#                     while a == 0:
#                         eigenValues, eigenVectors = np.linalg.eig(newMatrix)
#                         minEigen = np.min(abs(eigenValues))
#                         reg_ = min(abs(minEigen), 1e-6)
#                         newDiag = np.diag(eigenValues + reg_)
#                         newMatrix = eigenVectors @ newDiag @ eigenVectors.T
#                         try: 
#                             np.linalg.cholesky(newMatrix)
#                             a = 1
#                         except:
#                             a = 0
#                     with open(pathAnalysis + '/' + nameMatrix + '_correction.pkl','wb') as wf_:
#                         pickle.dump(newMatrix, wf_)
#                 else:
#                     with open(pathAnalysis + '/' + nameMatrix + '_correction.pkl','wb') as wf_:
#                         pickle.dump(desiredMatrix, wf_)
#             else:
#                 print(nameMatrix + ' already corrected!')
#     return 1

# Function equivalent to correct_K_C_alpha; Only difference is the name of input matrix
# def correct_GCTA(matrixes: dict, selectedSample: dict, check_GCTA_calculated:int):
#     pathAnalysis = selectedSample['pathAnalysis']
#     for nameMatrixIt_ in matrixes:
#         nameMatrix = 'GCTA_' + nameMatrixIt_
#         # Reads and stores as a variable the calculated matrix
#         if not exists(pathAnalysis + '/' + nameMatrix + '_correction.pkl'):
#             with open(pathAnalysis+'/' + nameMatrix + '.pkl','rb') as rf_:
#                 desiredMatrix = pickle.load(rf_)
#             try:
#                 choleskyDecomp_ = np.linalg.cholesky(desiredMatrix)
#                 a = 1
#             except:
#                 a = 0
#             if a == 0:
#                 newMatrix = desiredMatrix.copy()
#                 while a == 0:
#                     eigenValues, eigenVectors = np.linalg.eig(newMatrix)
#                     minEigen = np.min(abs(eigenValues))
#                     reg_ = abs(minEigen) + 1e-6
#                     newDiag = np.diag(eigenValues + reg_)
#                     newMatrix = eigenVectors @ newDiag @ eigenVectors.T
#                     try: 
#                         np.linalg.cholesky(newMatrix)
#                         a = 1
#                     except:
#                         a = 0
#                 with open(pathAnalysis + '/' + nameMatrix + '_correction.pkl','wb') as wf_:
#                     pickle.dump(newMatrix, wf_)
#             else:
#                 with open(pathAnalysis + '/' + nameMatrix + '_correction.pkl','wb') as wf_:
#                     pickle.dump(desiredMatrix, wf_)
#         else:
#             print(nameMatrix + ' already corrected!')
#     return 1

def readBinaryGCTA(file_):
    t_ = pd.read_csv(f'{file_}.grm.id', sep = '\t', header = None)
    n_ = t_.shape[0]
    biteSize = 4
    totalElements = int(n_ * (n_ + 1) / 2)
    values_ = []
    with open(f'{file_}.grm.bin','rb') as binFile_:
        for i_ in range(totalElements):
            values_.append(struct.unpack('f',binFile_.read(biteSize))[0])
    dfIndex = pd.DataFrame([x for x in range(n_)], columns = ['index_'])
    dfIndex.loc[:,'index_'] = dfIndex.index_ + 1
    dfIndex.loc[:, 'indexC_'] = dfIndex.index_.cumsum()
    dfIndex.loc[:, 'indexDiagonal'] = dfIndex.indexC_
    indexDiagonal = dfIndex.indexDiagonal.values - 1
    diagonalMatrix = [values_[x] for x in indexDiagonal]
    offDiagonalMatrix = [values_[x] for x in range(len(values_)) if x not in indexDiagonal]
    # Build matrix from binary files
    gctaMatrix_ = np.zeros([n_,n_])
    cont_ = 0
    contDiag_ = 0
    for col_ in range(n_):
        for row_ in range(n_):
            if row_ < col_:
                gctaMatrix_[row_,col_] = offDiagonalMatrix[cont_]
                gctaMatrix_[col_,row_] = offDiagonalMatrix[cont_]
                cont_ += 1
            if row_ == col_:
                gctaMatrix_[row_,col_] = diagonalMatrix[contDiag_]
                contDiag_ += 1
    return gctaMatrix_

def calculate_K_matrix(params_, selectedSample):
    pathAnalysis = selectedSample['pathAnalysis']
    matrices_params_ = params_['matrices']
    sets_ = matrices_params_['sets_']
    print(sets_)
    type_ = matrices_params_['type_']
    for case_ in sets_:
        print(case_)
        for list_ in sets_[case_]:
            print(sets_[case_], list_)
            chromosomes = sets_[case_][list_]
            for matrix_ in type_:
                logging.info(matrix_)
                if 'K_C' in matrix_:
                    for alpha_ in type_[matrix_]:
                        nameMatrix = f'''K_C_{alpha_}_{case_}_{list_}'''
                        if not exists(f'''{pathAnalysis}/{nameMatrix}.pkl'''):
                            params = {}
                            params["pathAnalysis"] = pathAnalysis
                            params["alpha"] = alpha_
                            params["chromosomes"] = chromosomes
                            calculator = calculateMatrixes(params)
                            calculator.readBeds()
                            calculator.calculateMatrix()
                            desiredMatrix = calculator.zztFinal
                            with open(f'''{pathAnalysis}/{nameMatrix}.pkl''','wb') as wf_:
                                pickle.dump(desiredMatrix, wf_)
                            df_desired_matrix = pd.DataFrame(desiredMatrix)
                            df_desired_matrix.to_csv(f'''{pathAnalysis}/{nameMatrix}.txt''', sep = '|', header = None)
                        else:
                            print(f"""{nameMatrix} already calculated!""")
                elif 'GCTA' in matrix_:
                    nameMatrix = f'''GCTA_{case_}_{list_}'''
                    if not exists(f'''{pathAnalysis}/{nameMatrix}.pkl'''):
                        # Specific parameters to build matrices
                        params = {}
                        params["pathAnalysis"] = pathAnalysis
                        params["chromosomes"] = chromosomes
                        # Number of threads
                        numThreads_ = type_[matrix_]
                        # Hiden parameters, to run GCTA software
                        conf_paths = ["conf/local"]
                        conf_loader = ConfigLoader(conf_paths)
                        parameters = conf_loader.get("paths*", "paths*/**")
                        pathGCTA_ = parameters['gcta']
                        # Write ref file
                        ref_file_name_ = f'{pathAnalysis}/REF_BEDS_{case_}_{list_}'
                        file_ = open(ref_file_name_, 'w')
                        [file_.write(f'{pathAnalysis}/pruned_chr{chr_}\n') for chr_ in chromosomes]
                        file_.close()
                        # command to build matrix
                        cmd = [pathGCTA_, '--mbfile' ,ref_file_name_, '--keep', f'{pathAnalysis}/sample.txt' ,'--make-grm','--out',f'{pathAnalysis}/{nameMatrix}','--thread-num',str(numThreads_)]
                        # Execute command and wait for the end
                        p = subprocess.Popen(cmd)
                        p.wait()
                        # Read matrix
                        desiredMatrix = readBinaryGCTA(f'{pathAnalysis}/{nameMatrix}')
                        with open(f'''{pathAnalysis}/{nameMatrix}.pkl''','wb') as wf_:
                            pickle.dump(desiredMatrix, wf_)
                        df_desired_matrix = pd.DataFrame(desiredMatrix)
                        df_desired_matrix.to_csv(f'''{pathAnalysis}/{nameMatrix}.txt''', sep = '|', header = None)
                    else:
                        print(f"""{nameMatrix} already calculated!""")
    return selectedSample

def standardize_matrix(matrix_):
    # Extract main diagonal from matrix and build diagonal matrix with the inverse of its sqrt
    diag_matrix = np.diag(1/np.sqrt(np.diag(matrix_)))
    std_matrix = diag_matrix @ matrix_ @ diag_matrix
    return std_matrix

def standardize_matrix_if_asked(params_, selectedSample):
    pathAnalysis = selectedSample['pathAnalysis']
    matrices_params_ = params_['matrices']
    sets_ = matrices_params_['sets_']
    type_ = matrices_params_['type_']
    for case_ in sets_:
        for list_ in sets_[case_]:
            for matrix_ in type_:
                logging.info(matrix_)
                if 'std_' in matrix_:
                    if 'K_C' in matrix_:
                        for alpha_ in type_[matrix_]:
                            nameMatrixOrig = f'''K_C_{alpha_}_{case_}_{list_}'''
                            nameMatrix = f'''{matrix_}_{alpha_}_{case_}_{list_}'''
                            if not exists(f'''{pathAnalysis}/{nameMatrix}.pkl'''):
                                original_matrix = pickle.load(open(f'''{pathAnalysis}/{nameMatrixOrig}.pkl''', 'rb'))
                                std_matrix = standardize_matrix(original_matrix)
                                with open(f'''{pathAnalysis}/{nameMatrix}.pkl''','wb') as wf_:
                                    pickle.dump(std_matrix, wf_)
                                df_desired_matrix = pd.DataFrame(std_matrix)
                                df_desired_matrix.to_csv(f'''{pathAnalysis}/{nameMatrix}.txt''', sep = '|', header = None)
                            else:
                                print(f"""{nameMatrix} already standardized!""")
                    elif 'GCTA' in matrix_:
                        nameMatrixOrig = f'''GCTA_{case_}_{list_}'''
                        nameMatrix = f'''{matrix_}_{case_}_{list_}'''
                        if not exists(f'''{pathAnalysis}/{nameMatrix}.pkl'''):
                            original_matrix = pickle.load(open(f'''{pathAnalysis}/{nameMatrixOrig}.pkl''', 'rb'))
                            std_matrix = standardize_matrix(original_matrix)
                            with open(f'''{pathAnalysis}/{nameMatrix}.pkl''','wb') as wf_:
                                pickle.dump(std_matrix, wf_)
                            df_desired_matrix = pd.DataFrame(std_matrix)
                            df_desired_matrix.to_csv(f'''{pathAnalysis}/{nameMatrix}.txt''', sep = '|', header = None)
                        else:
                            print(f"""{nameMatrix} already standardized!""")
    return selectedSample

def make_non_negative(input_matrix_):
    try:
        np.linalg.cholesky(input_matrix_)
        a = 1
    except:
        a = 0
    newMatrix = input_matrix_.copy()
    # If the matrix is not positive definite, then proceed to algorithm
    # Else, just skip that part and return the input
    if a == 0:
        logging.info(f"""
        The matrix is not positive semidefinite
        """)
        counter_ = 0
        while a == 0:
            counter_ += 1
            eigenValues, eigenVectors = np.linalg.eig(newMatrix)
            minEigen = np.min(abs(eigenValues)) + 1e-6
            reg_ = min(abs(minEigen), 1e-4)
            newDiag = np.diag(eigenValues + reg_)
            logging.info(f"""
            Increasing {reg_} to main diagonal
            Current minimum eigen value is: {minEigen}
            New minimum value is {min(eigenValues + reg_)}
            """)
            newMatrix = eigenVectors @ newDiag @ eigenVectors.T
            try: 
                np.linalg.cholesky(newMatrix)
                a = 1
            except:
                a = 0
    else:
        counter_ = 0
        logging.info(f"""
        The matrix was already positive semidefinite
        """)
    return newMatrix, counter_

def make_K_matrix_non_negative(params_, selectedSample):
    pathAnalysis = selectedSample['pathAnalysis']
    matrices_params_ = params_['matrices']
    sets_ = matrices_params_['sets_']
    type_ = matrices_params_['type_']
    for case_ in sets_:
        for list_ in sets_[case_]:
            for matrix_ in type_:
                if 'K_C' in matrix_:
                    for alpha_ in type_[matrix_]:
                        nameMatrix = f'''{matrix_}_{alpha_}_{case_}_{list_}'''
                        if not exists(f'''{pathAnalysis}/{nameMatrix}_correction_non_negative.txt'''):
                            with open(f'''{pathAnalysis}/{nameMatrix}.pkl''','rb') as rf_:
                                desiredMatrix = pickle.load(rf_)
                            newMatrix, counter_ = make_non_negative(desiredMatrix)
                            with open(f'''{pathAnalysis}/{nameMatrix}_correction_non_negative.pkl''','wb') as wf_:
                                pickle.dump(newMatrix, wf_)
                            newMatrix = pd.DataFrame(newMatrix)
                            newMatrix.to_csv(f'''{pathAnalysis}/{nameMatrix}_correction_non_negative.txt''', sep = '|', header = None)
                            logging.info(f"""
                            Done! Matrix is positive semidefinite
                            Number of trials until matrix {nameMatrix} correction: {counter_}
                            """)
                        else:
                            logging.info(f'''{nameMatrix} already corrected!''')
                elif 'GCTA' in matrix_:
                    nameMatrix = f'''{matrix_}_{case_}_{list_}'''
                    if not exists(f'''{pathAnalysis}/{nameMatrix}_correction_non_negative.txt'''):
                        desiredMatrix = pickle.load(open(f'''{pathAnalysis}/{nameMatrix}.pkl''', 'rb'))
                        newMatrix, counter_ = make_non_negative(desiredMatrix)
                        with open(f'''{pathAnalysis}/{nameMatrix}_correction_non_negative.pkl''','wb') as wf_:
                            pickle.dump(newMatrix, wf_)
                        newMatrix = pd.DataFrame(newMatrix)
                        newMatrix.to_csv(f'''{pathAnalysis}/{nameMatrix}_correction_non_negative.txt''', sep = '|', header = None)
                        logging.info(f"""
                        Done! Matrix is positive semidefinite
                        Number of trials until matrix {nameMatrix} correction: {counter_}
                        """)
                    else:
                        logging.info(f'''{nameMatrix} already corrected!''')
    return selectedSample



