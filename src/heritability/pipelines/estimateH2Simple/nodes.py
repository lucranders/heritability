"""
This is a boilerplate pipeline 'estimateH2Simple'
generated using Kedro 0.17.6
"""
import pandas as pd
from . import estimateH2 as herit
from kedro.config import ConfigLoader

def assembleParams(snpsParams: dict , sampParams: dict , formula: dict, GeneExpressions: list , pathTemp:str,  sampSelection: dict):
    args_ = dict()
    conf_paths = ["conf/local"]
    conf_loader = ConfigLoader(conf_paths)
    parameters = conf_loader.get("paths*", "paths*/**")
    saveFile = parameters['pathSaveFilesSimple'] 
    args_['vif_'] = snpsParams['vif']
    args_['maf_'] = snpsParams['maf']
    args_['hwe_'] = snpsParams['hwe']
    pop_ = ''
    for i,x in enumerate( sampParams['pop'] ):
        if i < ( len(sampParams['pop']) - 1 ):
            pop_ += x + '_'
        else:
            pop_ += x
    args_['pop_'] = pop_
    args_['outLiers_'] = sampParams['outliers']
    args_['sizeOriginalDf'] = sampSelection['sizeOriginalDf']
    args_['sizeFinalDf'] = sampSelection['sizeFinalDf']
    args_['formulaFE_'] = formula['fixed']
    args_['genes'] = GeneExpressions
    args_['tempFolder'] = pathTemp
    args_['db'] = sampSelection['originalDf']
    args_['individuals'] = sampSelection['selectedSample']
    args_['extraRE'] = formula['random']
    args_['saveRef_'] = saveFile
    return  args_

def estimateSigmas2REMLSimple(snpsParams: dict , sampParams: dict , formula: dict, GeneExpressions: list , pathTemp:str,  sampSelection: dict , genesMatrix: dict):
    generalParams_ = assembleParams(snpsParams , sampParams , formula, GeneExpressions , pathTemp,  sampSelection)
    print(generalParams_)
    generalParams_['additiveMatrixDictionary'] = {'Genes': genesMatrix}
    print(generalParams_['additiveMatrixDictionary'])
    generalParams_['method'] = 'REML'
    generalParams_['parametersOpt'] = {'maxIt':100,'conv':1e-4}
    for sigma0 in range(10000 , 100000 , 10000):
        sigmasInit = {'Residuals': sigma0 , 'Genes': 100000 - sigma0}
        generalParams_['parametersOpt']['sigmasInit'] = sigmasInit
        estimateH2_REMLSimple = herit.heritabilityAlt(generalParams_)
        estimateH2_REMLSimple.calculateAll('Genes')
        estimateH2_REMLSimple.saveResults()