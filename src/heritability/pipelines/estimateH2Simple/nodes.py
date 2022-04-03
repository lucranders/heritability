"""
This is a boilerplate pipeline 'estimateH2Simple'
generated using Kedro 0.17.6
"""
import pandas as pd
import numpy as np
import patsy

def defineFEM_PV(data: pd.DataFrame, geneExpr: str, pathTempFiles:str, formulaFE:str) -> pd.DataFrame:
    """
    This function is designed to build the Design Matrix (Fixed Effects) and the response variable
    (gene expression). For each gene expression, this function need to be executed
    """
    sampleToCalculateMatrix = pd.read_csv(pathTempFiles + '/sample.txt', header= None, sep = ' ')
    sampleToCalculateMatrix.columns = ['subject_id','drop']
    sampleToCalculateMatrix.drop(columns=['drop'], inplace = True)
    filterGene = data.loc[data.gene_name == geneExpr].copy()
    givenDb = sampleToCalculateMatrix.merge(filterGene)
    XPop_ = patsy.dmatrix( formulaFE ,data = givenDb)
    YPop_ = givenDb.loc[:,geneExpr].copy()
    return dict(xMatrix=XPop_, yVector=YPop_.values, filteredDf=givenDb)

# Define Random Effects Matrix (Covariances)
def defineREM(data:pd.DataFrame, additiveMatrixDictionary:dict, extraRE:list) -> dict:
    """
    This function is designed to build extra random effects Matrices, if desired
    If not, just adds the residuals matrix (Identity) to additiveMatrixDictionary
    dictionary. Parameters:
    
    data: dataframe from defineFEM_PV execution (filteredDf)
    additiveMatrixDictionary: Dictionary with previously calculated matrices (from previous
    steps)
    extraRE: List containing desired random effects
    """
    if extraRE != None:
        for var_ in extraRE:
            oneHotEncoding = pd.get_dummies(data.loc[:, var_]).values
            covMat = np.matmul( oneHotEncoding , np.transpose(oneHotEncoding) )
            additiveMatrixDictionary[var_] = covMat.copy()
            del(covMat)
    # Diagonal matrix, representing residuals
    additiveMatrixDictionary['Residuals'] = np.diag(np.full(additiveMatrixDictionary[list(additiveMatrixDictionary)[0]].shape[0],1))
# Estimate parameters from additive model through Maximum Likelihood
def estimateAddVarML(xMatrix: np.matrix, yVector:np.matrix, randomCovs:dict,parametersOpt):
    # main definitions to start algorithm
    X = xMatrix.copy()
    y = yVector.copy()
    givenMatrixes = randomCovs.copy()
    dimSquareMatrixes = givenMatrixes[list(givenMatrixes)[0]].shape[0]
    nameComponents = list(givenMatrixes)
    numberMatrixes = len(givenMatrixes)
    if 'sigmasInit' in list(parametersOpt):
        sigmas2Init = parametersOpt['sigmasInit']
    else:
        # if not given initial values, initializes with variance from expression, divided by number of covariance matrixes
        sigmas2Init = dict()
        for var_ in nameComponents: 
            sigmas2Init[var_] = y.var()/numberMatrixes
    if 'conv' in list(parametersOpt):
        conv = parametersOpt['conv']
    else:
        conv = 1e-4
    if 'maxIt' in list(parametersOpt):
        maxIt = parametersOpt['maxIt']
    else:
        maxIt = 10
    betas = dict()
    logLik = dict()
    # Initializing - first iteraction
    V = np.zeros([dimSquareMatrixes,dimSquareMatrixes])
    for var_ in nameComponents:
        V = V + (givenMatrixes[var_] * sigmas2Init[var_])
    vInverse = np.linalg.inv(V)
    beta_ = np.linalg.inv(np.transpose(X) @ vInverse @ X) @ np.transpose(X) @ vInverse @ y
    difMean = y - (X @ beta_)
    likelihood_ = np.linalg.slogdet(V)[1] + (np.transpose(difMean) @ vInverse @ difMean)
    logLik[0] = -.5*likelihood_.copy()
    betas[0] = beta_.copy()
    thetas = dict()
    thetas[0] = sigmas2Init.copy()
    k = 1
    while k < maxIt:
        s = list()
        Hessian = np.zeros([numberMatrixes,numberMatrixes])
        for numRow,cov1 in enumerate(nameComponents):
            # score vector element
            sElement = np.trace(vInverse @ givenMatrixes[cov1]) - (np.transpose(difMean) @ vInverse @ givenMatrixes[cov1] @ vInverse @ difMean)
            s.append(-.5*sElement.copy())
            del(sElement)
            for numCol,cov2 in enumerate(nameComponents):
                hElement = np.trace(vInverse @ givenMatrixes[cov1] @ vInverse @ givenMatrixes[cov2])
                Hessian[numRow,numCol] = -.5*hElement.copy()
                del(hElement)
        invHessian = np.linalg.inv(Hessian)
        V = np.zeros([dimSquareMatrixes,dimSquareMatrixes])
        thetas_ = dict()
        for numRow,var_ in enumerate(nameComponents):
            thetaNew = sigmas2Init[var_].copy() - (invHessian[numRow,:] @ np.c_[s])
            thetas_[var_] = thetaNew.copy()
            V = V + (givenMatrixes[var_] * thetas_[var_])
        sigmas2Init = thetas_.copy()
        thetas[k] = thetas_.copy()
        vInverse = np.linalg.inv(V)
        beta_ = np.linalg.inv(np.transpose(X) @ vInverse @ X) @ np.transpose(X) @ vInverse @ y
        betas[k] = beta_.copy()
        difMean = y - (X @ beta_)
        likelihood_ = np.linalg.slogdet(V)[1] + (np.transpose(difMean) @ vInverse @ difMean)
        logLik[k] = -.5*likelihood_.copy()
        print(logLik[k])
        if np.absolute(logLik[k] - logLik[k-1]) < conv:
            break
        k = k+1
    return dict(sigma2 = thetas, logLikelihood = logLik, betas = betas)
# Estimate parameters from additive model through Restricted/Residual Maximum Likelihood (REML)
def estimateAddVarREML(xMatrix: np.matrix, yVector:np.matrix, randomCovs:dict,parametersOpt):
        # main definitions to start algorithm
        X = xMatrix.copy()
        y = yVector.copy()
        givenMatrixes = randomCovs.copy()
        dimSquareMatrixes = givenMatrixes[list(givenMatrixes)[0]].shape[0]
        nameComponents = list(givenMatrixes)
        numberMatrixes = len(givenMatrixes)
        if 'sigmasInit' in list(parametersOpt):
            sigmas2Init = parametersOpt['sigmasInit']
        else:
            # if not given initial values, initializes with variance from expression, divided by number of covariance matrixes
            sigmas2Init = dict()
            for var_ in nameComponents: 
                sigmas2Init[var_] = y.var()/numberMatrixes
        if 'conv' in list(parametersOpt):
            conv = parametersOpt['conv']
        else:
            conv = 1e-4
        if 'maxIt' in list(parametersOpt):
            maxIt = parametersOpt['maxIt']
        else:
            maxIt = 10
        betas = dict()
        logLik = dict()
        # Initializing - first iteraction
        V = np.zeros([dimSquareMatrixes,dimSquareMatrixes])
        for var_ in nameComponents:
            V = V + (givenMatrixes[var_] * sigmas2Init[var_])
        vInverse = np.linalg.inv(V)
        beta_ = np.linalg.inv(np.transpose(X) @ vInverse @ X) @ np.transpose(X) @ vInverse @ y
        difMean = y - (X @ beta_)
        P = vInverse - ( vInverse @ X @ np.linalg.inv( np.transpose(X) @ vInverse @ X)  @ np.transpose(X) @ vInverse )
        likelihood_ = np.linalg.slogdet(V)[1] + np.linalg.slogdet( np.transpose(X) @ vInverse @ X )[1] + (np.transpose(difMean) @ vInverse @ difMean)
        logLik[0] = -.5*likelihood_.copy()
        betas[0] = beta_.copy()
        thetas = dict()
        thetas[0] = sigmas2Init.copy()
        k = 1
        while k < maxIt:
            s = list()
            Hessian = np.zeros([numberMatrixes,numberMatrixes])
            for numRow,cov1 in enumerate(nameComponents):
                # score vector element
                sElement = np.trace(P @ givenMatrixes[cov1]) - np.transpose(difMean) @ vInverse @ givenMatrixes[cov1] @ vInverse @ difMean
                s.append(-.5*sElement.copy())
                del(sElement)
                for numCol,cov2 in enumerate(nameComponents):
                    hElement = np.trace(P @ givenMatrixes[cov1] @ P @ givenMatrixes[cov2])
                    Hessian[numRow,numCol] = -.5*hElement.copy()
                    del(hElement)
            invHessian = np.linalg.inv(Hessian)
            V = np.zeros([dimSquareMatrixes,dimSquareMatrixes])
            thetas_ = dict()
            for numRow,var_ in enumerate(nameComponents):
                thetaNew = sigmas2Init[var_].copy() - (invHessian[numRow,:] @ np.c_[s])
                thetas_[var_] = thetaNew.copy()
                V = V + (givenMatrixes[var_] * thetas_[var_])
            vInverse = np.linalg.inv(V)
            P = vInverse - ( vInverse @ X @ np.linalg.inv( np.transpose(X) @ vInverse @ X)  @ np.transpose(X) @ vInverse )
            sigmas2Init = thetas_.copy()
            thetas[k] = thetas_.copy()
            beta_ = np.linalg.inv(np.transpose(X) @ vInverse @ X) @ np.transpose(X) @ vInverse @ y
            betas[k] = beta_.copy()
            difMean = y - (X @ beta_)
            likelihood_ = np.linalg.slogdet(V)[1] + np.linalg.slogdet( np.transpose(X) @ vInverse @ X )[1] + (np.transpose(difMean) @ vInverse @ difMean)
            logLik[k] = -.5*likelihood_.copy()
            print(logLik[k])
            if np.absolute(logLik[k] - logLik[k-1]) < conv:
                break
            k = k+1
        return dict(sigma2 = thetas, logLikelihood = logLik,betas = betas)
