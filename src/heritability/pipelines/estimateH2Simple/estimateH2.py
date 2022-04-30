import subprocess
import os
import time
from datetime import datetime
from pathlib import Path
import pandas as pd
import numpy as np
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
import patsy
from multiprocessing import Pool
import logging

def extractGCTAResults(path_):
    f_ = open(path_,'r')
    for line in f_:
        strings_ = line.split('\t')
        if 'Vp' in line and 'V(G)' not in line:
            totalVariance = strings_[1]
        if 'Vp' not in line and 'V(G)' in line:
            varianceSnps = strings_[1]
        if 'Vp' in line and 'V(G)' in line:
            herit = strings_[1]
    return float(varianceSnps) , float(totalVariance) , float(herit)

class heritabilityGCTA:
    def __init__(self,individuals,db,covs,genes,oldParams):
        self.pop_ = oldParams.pop_
        self.maf_ = oldParams.maf_
        self.hwe_ = oldParams.hwe_
        self.vif_ = oldParams.vif_
        self.genes_ = genes
        self.covs_ = covs
        self.pathPipeline_ = oldParams.pathPipeline_
        self.pathTmp_ = oldParams.pathTmp_
        self.pathGCTA_ = oldParams.pathGCTA_
        self.path_ = self.pathTmp_ + '/TempBed_pop_' + self.pop_ + "_maf_" + str(self.maf_) + "_hwe_" + str(self.hwe_) + "_vif_" + str(self.vif_)
        self.db_ = db
        self.individuals_ = individuals
    # Creates required files to estimate heritability through GCTA software
    def createFiles(self,geneExpr_):
        basicCols = ['subject_id']
        X_ = basicCols.copy()
        Y_ = basicCols.copy()
        Y_.append('tpm')
        print(self.db_.head(3))
        for cov_ in self.covs_:
            X_.append(cov_)
        dbToMerge = self.individuals_.copy()
        covariatesX = self.db_.loc[self.db_.gene_name == geneExpr_,X_]
        expressionY = self.db_.loc[self.db_.gene_name == geneExpr_,Y_]
        XPop_ = pd.concat([dbToMerge.loc[:,['subject_id']] , dbToMerge.merge(covariatesX)] , axis = 1)
        YPop_ = pd.concat([dbToMerge.loc[:,['subject_id']] , dbToMerge.merge(expressionY)] , axis = 1)
        XPop_.to_csv(self.path_ + '/X_' + geneExpr_.replace('-','') +'.txt',sep = ' ',index = False, header = False)
        YPop_.to_csv(self.path_ + '/Y_' + geneExpr_.replace('-','') +'.txt',sep = ' ',index = False, header = False)
    # Estimates heritability through GCTA software for one gene expression
    def heritabilityGCTA(self,geneExpr_,suffix_ = None):
        gctaLoc = '/raid/genevol/users/lucas/gcta/gcta64'
        method = '--reml'
        algorithm = '--reml-alg'
        algorithmV = '0'
        maxit = '--reml-maxit'
        maxitV = '1000'
        grm = '--grm'
        if suffix_ == None:
            grmV = self.path_ + '/GCTA'
        else:
            grmV = self.path_ + '/GCTA' + suffix_
        pheno =  '--pheno'
        phenoV = self.path_ + '/Y_' + geneExpr_.replace('-','') + '.txt'
        covs = '--covar'
        covsV = self.path_ + '/X_' + geneExpr_.replace('-','') + '.txt'
        fix = '--reml-est-fix'
        out = '--out'
        outV = self.path_ + '/Results' + geneExpr_.replace('-','')
        cmd = [gctaLoc, method, algorithm, algorithmV, maxit, maxitV, grm, grmV, pheno, phenoV, covs, covsV, fix, out, outV]
        subprocess.Popen(cmd)
    # Estimates heritability through GCTA software for all gene expressions
    def calculateAll(self,suffix_):
        for geneExpr_ in self.genes_:
            self.createFiles(geneExpr_)
            self.heritabilityGCTA(geneExpr_,suffix_)
    def resultsToDf(self,sampInf_,formulaFE_,formulaRE_ = 'Genes+Residuals'):
        finalTuple = []
        for geneExpr_ in self.genes_:
            gene_ = geneExpr_.replace('-','')
            pathGenes_ = self.path_ + '/Results' + gene_ + '.log'
            snpsGCTA, totalGCTA, h2GCTA = extractGCTAResults(pathGenes_)
            pop = self.pop_
            maf = self.maf_
            hwe = self.hwe_
            vif = self.vif_
            sampInf = sampInf_
            tuple_ = [geneExpr_,pop , maf , hwe , vif , sampInf,h2GCTA,snpsGCTA,totalGCTA,'GCTA',formulaFE_,formulaRE_]
            finalTuple.append(tuple_)
        finalDf = pd.DataFrame(finalTuple)
        finalDf.columns = ['Gene','Pop' , 'MAF' , 'HWE' , 'VIF' , 'sampleInference','HeritEst','DesiredVarEst','totalVarEst','method','formulaF','randEffects']
        self.results = finalDf
    def saveResults(self,saveRef_):
        check_ = 0
        try:
            df_ = pd.read_csv(saveRef_,sep = '|',header = 0)
        except:
            check_ = 1
        if check_ == 1:
            self.results.to_csv(saveRef_,index = False, header = True, sep = "|")
        else:
            df_ = pd.concat([df_,self.results], axis = 0)
            df_.drop_duplicates(inplace = True)
            df_.to_csv(saveRef_,index = False, header = True, sep = "|")


class heritabilityAlt:
    def __init__(self,args):
        self.pop_ = args['pop_']
        self.maf_ = args['maf_']
        self.hwe_ = args['hwe_']
        self.vif_ = args['vif_']
        self.outLiers_ = args['outLiers_']
        self.sizeOriginalDf = args['sizeOriginalDf']
        self.sizeFinalDf = args['sizeFinalDf']
        self.formulaFE_ = args['formulaFE_']
        self.genes_ = args['genes']
        self.path_ = args['tempFolder']
        self.db_ = args['db']
        self.individuals_ = args['individuals']
        self.additiveMatrixDictionary_ = args['additiveMatrixDictionary']
        self.extraRE_ = args['extraRE']
        self.parametersOpt_ = args['parametersOpt']
        self.method_ = args['method']
        self.saveRef_ = args['saveRef_']
    # Define Fixed Effects Matrix and Phenotype Vector
    def defineFEM_PV(self,geneExpr_):
        dbToMerge = self.individuals_.copy()
        filterGene = self.db_.loc[self.db_.gene_name == geneExpr_].copy()
        givenDb = dbToMerge.merge(filterGene)
        XPop_ = patsy.dmatrix( self.formulaFE_ ,data = givenDb)
        YPop_ = givenDb.loc[:,['gene_name']].copy()
        self.xMatrix = XPop_
        self.yVector = YPop_.values
    # Define Random Effects Matrix (Covariances)
    def defineREM(self,geneExpr_):
        self.randomCovs = self.additiveMatrixDictionary_
        if self.extraRE_ != None:
            dbToMerge = self.individuals_.copy()
            filterGene = self.db_.loc[self.db_.gene_name == geneExpr_].copy()
            givenDb = dbToMerge.merge(filterGene)
            for var_ in self.extraRE_:
                print(var_)
                oneHotEncoding = pd.get_dummies(givenDb.loc[:, var_]).values
                covMat = np.matmul( oneHotEncoding , np.transpose(oneHotEncoding) )
                self.randomCovs[var_] = covMat.copy()
                del(covMat)
        # Diagonal matrix, representing residuals
        self.randomCovs['Residuals'] = np.diag(np.full(self.additiveMatrixDictionary_[list(self.additiveMatrixDictionary_)[0]].shape[0],1))
    # Estimate parameters from additive model through Maximum Likelihood
    def estimateAddVarML(self):
        # main definitions to start algorithm
        X = self.xMatrix.copy()
        y = self.yVector.copy()
        givenMatrixes = self.randomCovs.copy()
        dimSquareMatrixes = givenMatrixes[list(givenMatrixes)[0]].shape[0]
        nameComponents = list(givenMatrixes)
        numberMatrixes = len(givenMatrixes)
        if 'sigmasInit' in list(self.parametersOpt_):
            sigmas2Init = self.parametersOpt_['sigmasInit']
        else:
            # if not given initial values, initializes with variance from expression, divided by number of covariance matrixes
            sigmas2Init = dict()
            for var_ in nameComponents: 
                sigmas2Init[var_] = y.var()/numberMatrixes
        if 'conv' in list(self.parametersOpt_):
            conv = self.parametersOpt_['conv']
        else:
            conv = 1e-4
        if 'maxIt' in list(self.parametersOpt_):
            maxIt = self.parametersOpt_['maxIt']
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
        self.sigma2 = thetas
        self.logLikelihood = logLik
        self.betas = betas
    # Estimate parameters from additive model through Restricted/Residual Maximum Likelihood (REML)
    def estimateAddVarREML(self):
        # main definitions to start algorithm
        X = self.xMatrix.copy()
        y = self.yVector.copy()
        givenMatrixes = self.randomCovs.copy()
        dimSquareMatrixes = givenMatrixes[list(givenMatrixes)[0]].shape[0]
        nameComponents = list(givenMatrixes)
        numberMatrixes = len(givenMatrixes)
        if 'sigmasInit' in list(self.parametersOpt_):
            sigmas2Init = self.parametersOpt_['sigmasInit']
        else:
            # if not given initial values, initializes with variance from expression, divided by number of covariance matrixes
            sigmas2Init = dict()
            for var_ in nameComponents: 
                sigmas2Init[var_] = y.var()/numberMatrixes
        if 'conv' in list(self.parametersOpt_):
            conv = self.parametersOpt_['conv']
        else:
            conv = 1e-4
        if 'maxIt' in list(self.parametersOpt_):
            maxIt = self.parametersOpt_['maxIt']
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
        self.sigma2 = thetas
        self.logLikelihood = logLik
        self.betas = betas
    def estimateSimpleHeritR(self,pathPipeline,geneExpr_,matrixName_,method_,fileSave_):
        cmd = ['Rscript',pathPipeline + '/calculateSimpleHeritR.r',self.path_,geneExpr_,matrixName_,method_,fileSave_]
        process = subprocess.Popen(cmd)
        process.wait()
    def calculateAllR(self,pathPipeline,matrixName_,method_,fileSave_):
        for geneExpr_ in self.genes_:
            self.estimateSimpleHeritR(pathPipeline = pathPipeline ,geneExpr_ = geneExpr_,matrixName_ = matrixName_,method_=method_,fileSave_=fileSave_)
    def resultsToDf(self,fileSave_,sampInf_,method_,formulaFE_,formulaRE_ = 'Genes+Residuals'):
        finalTuple = []
        for geneExpr_ in self.genes_:
            gene_ = geneExpr_.replace('-','')
            pathGenes_ = self.path_ + '/' + fileSave_ + '_' + gene_ + '.txt'
            snpsR, totalR, h2R = extractGCTAResults(pathGenes_)
            pop = self.pop_
            maf = self.maf_
            hwe = self.hwe_
            vif = self.vif_
            sampInf = sampInf_
            tuple_ = [geneExpr_,pop , maf , hwe , vif , sampInf,h2R,snpsR,totalR,'REstimate'+method_,formulaFE_,formulaRE_]
            finalTuple.append(tuple_)
        finalDf = pd.DataFrame(finalTuple)
        finalDf.columns = ['Gene','Pop' , 'MAF' , 'HWE' , 'VIF' , 'sampleInference','HeritEst','DesiredVarEst','totalVarEst','method','formulaF','randEffects']
        self.results = finalDf
    # Calculate heritability for one gene expression
    def getSigmas(self):
        totSum = 0
        lastRun = self.sigma2[len(self.sigma2)-1]
        for cov_ in list(lastRun):
            totSum += lastRun[cov_]
        self.finalDesiredSigma2 = lastRun
        self.finalTotalSigma = totSum
    # Calculate heritability for all gene expressions
    def calculateAll(self,column_):
        finalTuple = []
        for geneExpr_ in self.genes_:
            print(geneExpr_)
            self.defineFEM_PV(geneExpr_)
            self.defineREM(geneExpr_)
            if self.method_ == 'REML':
                try:
                    self.estimateAddVarREML()
                    a = 1
                except:
                    a = 0
            elif self.method_ == 'ML':
                try:
                    self.estimateAddVarML()
                    a = 1
                except:
                    a = 0
            else:
                print("Invalid option")
                break
            if a == 1:
                print('Success! Variances for ' + geneExpr_ + ' estimated')
                pop = self.pop_
                maf = self.maf_
                hwe = self.hwe_
                vif = self.vif_
                outliers = self.outLiers_
                sizeOriginalDf = self.sizeOriginalDf
                sizeFinalDf = self.sizeFinalDf
                namesRE = list(self.randomCovs)
                for num_,x in enumerate(namesRE):
                    if num_ == 0:
                        RE = x
                    else:
                        RE += "+" + x
                tuple_ = [geneExpr_,pop , maf , hwe , vif , outliers , sizeOriginalDf , sizeFinalDf ,self.method_,self.formulaFE_,RE]
                for sigmas in self.parametersOpt_:
                    tuple_.append(self.parametersOpt_[sigmas])
                self.getSigmas()
                finalDesiredSigma2 = self.finalDesiredSigma2[0].copy()
                finalTotalSigma = self.finalTotalSigma[0].copy()
                for cov_ in list(self.finalDesiredSigma2):
                    tuple_.append(self.finalDesiredSigma2[cov_][0].copy())
                finalTuple.append(tuple_)
                finalDf = pd.DataFrame(finalTuple)
                cols_ = ['Gene','Pop' , 'MAF' , 'HWE' , 'VIF' , 'outliers','sizeOriginalDf','sizeFinalDf','method','formulaF','randEffects']
                for sigmas in self.parametersOpt_:
                    cols_.append(sigmas+'_init')
                for cov_ in list(self.finalDesiredSigma2):
                    cols_.append(cov_+'_final')
                finalDf.columns = cols_
                self.results = finalDf
            else:
                print('Fail! Variances for ' + geneExpr_ + ' not estimated')
    # Save results
    def saveResults(self):
        check0_ = 0
        try:
            t_ = self.results
        except:
            check0_ = 1
        if check0_ == 0:
            check_ = 0
            try:
                df_ = pd.read_csv(self.saveRef_,sep = '|',header = 0)
            except:
                check_ = 1
            if check_ == 1:
                logging.info('A new file is created to store results: ' + self.saveRef_)
                self.results.to_csv(self.saveRef_,index = False, header = True, sep = "|")
            else:
                logging.info('The results will be appended to: ' + self.saveRef_)
                df_ = pd.concat([df_,self.results], axis = 0)
                df_.drop_duplicates(inplace = True)
                df_.to_csv(self.saveRef_,index = False, header = True, sep = "|")
        else:
            logging.error('No calculations were done!')
class selectNested:
    def __init__(self, dfPop, varPop, varSubject, maf, hwe, vif, pathTmp, pathPipeline):
        self.dfPop_ = dfPop
        self.varPop_ = varPop
        self.varSubject_ = varSubject
        self.difPops_ = self.dfPop_.loc[:,self.varPop_].unique()
        self.maf_ = maf
        self.vif_ = vif
        self.hwe_ = hwe
        self.path_ = dict()
        self.pathFilt_ = dict()
        self.pathPipeline_ = pathPipeline
        self.pathTmp_ = pathTmp
    # Create one temporary folder for each population
    def createTmpFolders(self):
        for pop_ in self.difPops_:
            self.path_[pop_] = self.pathTmp_ + '/TempBed_pop_' + pop_ + "_maf_" + str(self.maf_) + "_hwe_" + str(self.hwe_) + "_vif_" + str(self.vif_)
            subprocess.Popen(['mkdir',self.path_[pop_]])
            self.pathFilt_[pop_] = self.path_[pop_] + '/' + pop_ + '.filt'
            self.dfPop_.loc[self.dfPop_[self.varPop_] == pop_,[self.varSubject_,self.varSubject_]].drop_duplicates().to_csv(self.pathFilt_[pop_],index = False, header = False,sep = '\t')
    # Create list of snps from given parameters (maf, hwe and vif) and .bed files for each population, parallelized
    def createBedFileNPC(self):
        cmd = list()
        pool = Pool()
        for pop_ in self.difPops_:
            query_ = "plink --vcf $input"
            if self.vif_ != None:
                query_ += " --indep 50 5 " + str(self.vif_)
            if self.maf_ != None:
                query_ += " --maf " + str(self.maf_)
            if self.hwe_ != None:
                query_ += " --hwe " + str(self.hwe_)
            query_ += " --keep " + self.pathFilt_[pop_] + ' --mind 0.05 --geno 0.05 --make-bed --out $bed'
            # Create .bed files - 22 chromossomes
            cmd.append(['qsub', '-v' ,'query_=' + query_ + ',tempPath=' + self.path_[pop_] , self.pathPipeline_ + '/01.DataPrepGeneralNPC.sh'])
        pool.map(subprocess.Popen, cmd)
    # Writes files, containing status of Bed files creation for a given population
    def constantCheck(self,pop_):
        sizes = checkLogSizes(self.path_[pop_])
        cond = sum(sizes) < 22
        # If not, waits until so
        if cond == True:
            while cond:
                sizes = checkLogSizes(self.path_[pop_])
                cond = sum(sizes) < 22
                if cond:
                    # Updates status (inside tmp folder)
                    updateLog(self.path_[pop_],0,"bedStatus.txt")
                    time.sleep(10)
                else:
                    # Updates status (inside tmp folder)
                    updateLog(self.path_[pop_],1,"bedStatus.txt")
            # Create references to calculate GCTA
            for chr_ in range(1,23):
                f = open(self.path_[pop_] + '/chrs.txt', "a")
                f.write(self.path_[pop_]+'/chr' + str(chr_) + '\n')
                f.close()
    # Writes files, containing status of Bed files creation for all given populations, parallelized
    def checkBeds(self):
        pool = Pool()
        pool.map(self.constantCheck, list(self.difPops_))
