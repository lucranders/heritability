"""
This is a boilerplate pipeline 'estimateH2Simple'
generated using Kedro 0.17.6
"""
import pandas as pd
from . import estimateH2 as herit
from kedro.config import ConfigLoader
from os.path import exists
import pickle

def assembleParams(snpsParams: dict , sampParams: dict , formula: dict, GeneExpressions: list , selectedSample:dict, alphasExplored: dict, matrixes:dict):
    pathTemp = selectedSample['pathAnalysis']
    args_ = dict()
    conf_paths = ["conf/local"]
    conf_loader = ConfigLoader(conf_paths)
    parameters = conf_loader.get("paths*", "paths*/**")
    # There are different file layouts for each analysis (1, 2 or 22 sets of chromosomes)
    paramUse = {1:'pathSaveFilesSimple',2:'pathSaveFilesSimplePartitioned2Components',22:'pathSaveFilesSimplePartitioned22Components'}
    saveFile = parameters[paramUse[len(matrixes)]] 
    print(saveFile)
    args_['vif_'] = snpsParams['vif']
    args_['maf_'] = snpsParams['maf']
    args_['hwe_'] = snpsParams['hwe']
    pop_ = ''
    if sampParams['pop'] != None:
        for i,x in enumerate( sampParams['pop'] ):
            if i < ( len(sampParams['pop']) - 1 ):
                pop_ += x + '_'
            else:
                pop_ += x
    else:
        pop_ = 'All'
    args_['pop_'] = pop_
    args_['outLiers_'] = sampParams['outliers']
    args_['sizeOriginalDf'] = selectedSample['sizeOriginalDf']
    args_['sizeFinalDf'] = selectedSample['sizeFinalDf']
    args_['formulaFE_'] = formula['fixed']
    args_['genes'] = GeneExpressions
    args_['tempFolder'] = selectedSample['pathAnalysis']
    args_['db'] = selectedSample['originalDf']
    args_['individuals'] = selectedSample['selectedSample']
    args_['extraRE'] = formula['random']
    args_['saveRef_'] = saveFile
    args_['alphas'] = alphasExplored['alphas']
    return  args_

def estimateSigmas2REMLSimple(snpsParams: dict , sampParams: dict , formula: dict, GeneExpressions: list ,  selectedSample: dict , genesMatrixes: dict , saveControl: str):
    generalParams_ = assembleParams(snpsParams , sampParams , formula, GeneExpressions, selectedSample, genesMatrixes)
    # print(generalParams_['additiveMatrixDictionary'])
    generalParams_['method'] = 'REML'
    generalParams_['parametersOpt'] = {'maxIt':100,'conv':1e-4}
    if not exists(saveControl + '/done_'):
        for transf_ in ['None', 'log2', 'log10', 'ln']:
            for sigmaResid in range(10000 , 100001 , 10000):
                for sigmaGenes in range(10000 , 100001 , 10000):
                    generalParams_['transf_'] = transf_
                    sigmasInit = {'Residuals': sigmaResid , 'Genes': sigmaGenes}
                    generalParams_['parametersOpt']['sigmasInit'] = sigmasInit
                    estimateH2_REMLSimple = herit.heritabilityAlt(generalParams_)
                    estimateH2_REMLSimple.calculateAll('Genes')
                    estimateH2_REMLSimple.saveResults()
    else:
        print('calculations already done! To rerun, erase "done" file!')
    return 1

def estimateSigmas2REMLSimpleSingleStart(snpsParams: dict , sampParams: dict , formula: dict, GeneExpressions: list ,  selectedSample: dict , alphasExplored: dict, correctedMatrixes:dict, matrixes:dict):
    generalParams_ = assembleParams(snpsParams , sampParams , formula, GeneExpressions, selectedSample, alphasExplored, matrixes)
    generalParams_['method'] = 'REML'
    generalParams_['parametersOpt'] = {'maxIt':100,'conv':1e-4}
    pathAnalysis = selectedSample['pathAnalysis']
    for alpha_ in alphasExplored["alphas"]:
        dictGenes = {}
        for nameMatrixIt_ in matrixes:
            nameMatrix = 'K_C_' + str(alpha_) + '_' + nameMatrixIt_ 
            with open(pathAnalysis+'/' + nameMatrix + '_correction.pkl','rb') as rf_:
                desiredMatrix = pickle.load(rf_)
            dictGenes[nameMatrixIt_] = desiredMatrix
        generalParams_['additiveMatrixDictionary'] = dictGenes
        generalParams_['alpha'] = alpha_
        for transf_ in ['None', 'log2']:
            generalParams_['transf_'] = transf_
            estimateH2_REMLSimple = herit.heritabilityAlt(generalParams_)
            estimateH2_REMLSimple.calculateAll('Genes')
            estimateH2_REMLSimple.saveResults()
    return 1

def estimateSigmas2REMLSimpleSingleStartGCTA(snpsParams: dict , sampParams: dict , formula: dict, GeneExpressions: list ,  selectedSample: dict , alphasExplored: dict, correctedMatrixes:dict, matrixes:dict):
    generalParams_ = assembleParams(snpsParams , sampParams , formula, GeneExpressions, selectedSample, alphasExplored, matrixes)
    generalParams_['method'] = 'REML'
    generalParams_['parametersOpt'] = {'maxIt':100,'conv':1e-4}
    pathAnalysis = selectedSample['pathAnalysis']
    dictGenes = {}
    for nameMatrixIt_ in matrixes:
        nameMatrix = 'GCTA_' + nameMatrixIt_ 
        with open(pathAnalysis+'/' + nameMatrix + '_correction.pkl','rb') as rf_:
            desiredMatrix = pickle.load(rf_)
        dictGenes[nameMatrixIt_] = desiredMatrix
    generalParams_['additiveMatrixDictionary'] = dictGenes
    generalParams_['alpha'] = 'GCTA'
    for transf_ in ['None', 'log2']:
        generalParams_['transf_'] = transf_
        estimateH2_REMLSimple = herit.heritabilityAlt(generalParams_)
        estimateH2_REMLSimple.calculateAll('Genes')
        estimateH2_REMLSimple.saveResults()
    return 1

def execAllSimple(flgs1,flgs2,flgs3, saveControl: str):
    f = open(saveControl + '/done_','w')
    f.write('done!')
    f.close()
    return 1
# def estimateSigmas2MLSimple(snpsParams: dict , sampParams: dict , formula: dict, GeneExpressions: list , pathTemp:str,  sampSelection: dict , genesMatrix: dict):
#     generalParams_ = assembleParams(snpsParams , sampParams , formula, GeneExpressions , pathTemp,  sampSelection)
#     generalParams_['additiveMatrixDictionary'] = {'Genes': genesMatrix}
#     generalParams_['method'] = 'ML'
#     generalParams_['parametersOpt'] = {'maxIt':100,'conv':1e-4}
#     for sigmaResid in range(10000 , 100001 , 10000):
#         for sigmaGenes in range(10000 , 100001 , 10000):
#             sigmasInit = {'Residuals': sigmaResid , 'Genes': sigmaGenes}
#             generalParams_['parametersOpt']['sigmasInit'] = sigmasInit
#             estimateH2_REMLSimple = herit.heritabilityAlt(generalParams_)
#             estimateH2_REMLSimple.calculateAll('Genes')
#             estimateH2_REMLSimple.saveResults()
#     return 1
# def estimateSigmas2MLSimpleR(snpsParams: dict , sampParams: dict , formula: dict, GeneExpressions: list , pathTemp:str,  sampSelection: dict , genesMatrix: dict ):
#     generalParams_ = assembleParams(snpsParams , sampParams , formula, GeneExpressions , pathTemp,  sampSelection)
#     # print(generalParams_)
#     generalParams_['additiveMatrixDictionary'] = {'Genes': genesMatrix}
#     # print(generalParams_['additiveMatrixDictionary'])
#     generalParams_['method'] = 'ML'
#     generalParams_['parametersOpt'] = {'maxIt':100,'conv':1e-4}
#     for sigmaResid in range(10000 , 100001 , 10000):
#         for sigmaGenes in range(10000 , 100001 , 10000):
#             sigmasInit = {'Residuals': sigmaResid , 'Genes': sigmaGenes}
#             generalParams_['parametersOpt']['sigmasInit'] = sigmasInit
#             estimateH2_REMLSimple = herit.heritabilityAlt(generalParams_)
#             estimateH2_REMLSimple.calculateAllR('ML_','GCTA_Full_correction','testML')
#             estimateH2_REMLSimple.saveResults()
#     return 1