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


# Check log sizes from bash processes
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
class buildAdditiveVarianceMatrix:
    def __init__(self,samplesFile,pop,maf,hwe,vif,threads,pathPipeline,pathTmp,pathGCTA):
        self.pop_ = pop
        self.samplesFile_ = samplesFile
        self.maf_ = maf
        self.hwe_ = hwe
        self.vif_ = vif
        self.threads_ = str(threads)
        self.pathPipeline_ = pathPipeline
        self.pathTmp_ = pathTmp
        self.pathGCTA_ = pathGCTA
        self.path_ = pathTmp + '/TempBed_pop_' + pop + "_maf_" + str(maf) + "_hwe_" + str(hwe) + "_vif_" + str(vif)
    # Create .bed files
    def createTmpFolder(self):
        # Create temp folder
        subprocess.Popen(['mkdir',self.path_])
    def createBedFileNPC(self):
        query_ = "plink --vcf $input"
        if self.vif_ != None:
            query_ += " --indep 50 5 " + str(self.vif_)
        if self.maf_ != None:
            query_ += " --maf " + str(self.maf_)
        if self.hwe_ != None:
            query_ += " --hwe " + str(self.hwe_)
        if self.samplesFile_ != None:
            query_ += " --keep " + self.samplesFile_
        query_ += ' --mind 0.05 --geno 0.05 --make-bed --out $bed'
        self.query_ = query_
        # Create .bed files - 22 chromossomes
        cmd = ['qsub', '-v' ,'query_=' + self.query_ + ',tempPath=' + self.path_ , self.pathPipeline_ + '/01.DataPrepGeneralNPC.sh']
        subprocess.Popen(cmd)
        # Check if all 22 .bed files are created
        sizes = checkLogSizes(self.path_)
        cond = sum(sizes) < 22
        # If not, waits until so
        while cond:
            sizes = checkLogSizes(self.path_)
            cond = sum(sizes) < 22
            if cond:
                # Updates status (inside tmp folder)
                updateLog(self.path_,0,"bedStatus.txt")
                time.sleep(10)
            else:
                # Updates status (inside tmp folder)
                updateLog(self.path_,1,"bedStatus.txt")
        # Create references to calculate GCTA
        for chr_ in range(1,23):
            f = open(self.path_ + '/chrs.txt', "a")
            f.write(self.path_+'/chr' + str(chr_) + '\n')
            f.close()
    def createBedFile(self):
        # Create .bed files - 22 chromossomes
        cmd = ['qsub', '-v' ,'pop=' + self.pop_ + ',maf_=' + str(self.maf_) +',hwe_=' + str(self.hwe_) + ',vif_=' + str(self.vif_) + ',sampleFile_=' + str(self.sampleFile_) , self.pathPipeline_ + '/01.DataPrepGeneral.sh']
        subprocess.Popen(cmd)
        # Check if all 22 .bed files are created
        sizes = checkLogSizes(self.path_)
        cond = sum(sizes) < 22
        # If not, waits until so
        while cond:
            sizes = checkLogSizes(self.path_)
            cond = sum(sizes) < 22
            if cond:
                # Updates status (inside tmp folder)
                updateLog(self.path_,0,"bedStatus.txt")
                time.sleep(10)
            else:
                # Updates status (inside tmp folder)
                updateLog(self.path_,1,"bedStatus.txt")
        # Create references to calculate GCTA
        for chr_ in range(1,23):
            f = open(self.path_ + '/chrs.txt', "a")
            f.write(self.path_+'/chr' + str(chr_) + '\n')
            f.close()
    def calculateGCTA(self,refChrs,nameMatrix):
        if nameMatrix != None:
            self.nameMatrix_ = 'GCTA_' + nameMatrix
        else:
            self.nameMatrix_ = 'GCTA'
        if self.samplesFile_ != None:
            cmd = [self.pathGCTA_ + '/gcta64', '--mbfile' ,self.path_ + refChrs,'--keep',self.samplesFile_,'--make-grm','--out',self.path_+self.nameMatrix_,'--thread-num',self.threads_]
        else:
            cmd = [self.pathGCTA_ + '/gcta64', '--mbfile' ,self.path_ + refChrs,'--make-grm','--out',self.path_+self.nameMatrix_,'--thread-num',self.threads_]
        subprocess.Popen(cmd)
        # check whether GRM binaries already exists - means process is finished
        check_ = Path(self.path_ + self.nameMatrix_ + '.grm.bin')
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
    def correctGRM(self):
        subprocess.Popen(['Rscript',self.pathPipeline_ + '/02.matrixCorrection.r',self.path_,self.nameMatrix_])
        check_ = Path(self.path_ + '/statusGRMCorrection.txt')
        cond = check_.is_file()
        while not cond:
            cond = check_.is_file()
            time.sleep(1)
        subprocess.Popen(['rm',self.path_ + '/statusGRMCorrection.txt'])
    def readGRM(self):
        pandas2ri.activate()
        readRDS = robjects.r['readRDS']
        self.GRM = readRDS(self.path_ + '/' + self.nameMatrix_ +'.rds')
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
    def heritabilityGCTA(self,geneExpr_):
        gctaLoc = '/raid/genevol/users/lucas/gcta/gcta64'
        method = '--reml'
        algotithm = '--reml-alg'
        algotithmV = '0'
        maxit = '--reml-maxit'
        maxitV = '1000'
        grm = '--grm'
        grmV = self.path_ + '/GCTA'
        pheno =  '--pheno'
        phenoV = self.path_ + '/Y_' + geneExpr_.replace('-','') + '.txt'
        covs = '--covar'
        covsV = self.path_ + '/X_' + geneExpr_.replace('-','') + '.txt'
        fix = '--reml-est-fix'
        out = '--out'
        outV = self.path_ + '/Results' + geneExpr_.replace('-','')
        cmd = [gctaLoc, method, algotithm, algotithmV, maxit, maxitV, grm, grmV, pheno, phenoV, covs, covsV, fix, out, outV]
        subprocess.Popen(cmd)
    def calculateAll(self):
        for geneExpr_ in self.genes_:
            self.createFiles(geneExpr_)
            self.heritabilityGCTA(geneExpr_)
class heritabilityAlt:
    def __init__(self,individuals,sampInf,db,formulaFE,expression,genes,oldParams,additiveMatrixList,extraRE,parametersOpt,method,resultsFile):
        self.pop_ = oldParams.pop_
        self.maf_ = oldParams.maf_
        self.hwe_ = oldParams.hwe_
        self.vif_ = oldParams.vif_
        self.sampInf_ = sampInf
        self.formulaFE_ = formulaFE
        self.expression_ = expression
        self.genes_ = genes
        self.pathPipeline_ = oldParams.pathPipeline_
        self.pathTmp_ = oldParams.pathTmp_
        self.pathGCTA_ = oldParams.pathGCTA_
        self.path_ = self.pathTmp_ + '/TempBed_pop_' + self.pop_ + "_maf_" + str(self.maf_) + "_hwe_" + str(self.hwe_) + "_vif_" + str(self.vif_)
        self.db_ = db
        self.individuals_ = individuals
        self.additiveMatrixDictionary_ = additiveMatrixList
        self.extraRE_ = extraRE
        self.parametersOpt_ = parametersOpt
        self.method_ = method
        self.saveRef_ = resultsFile
    def defineFEM_PV(self,geneExpr_):
        dbToMerge = self.individuals_.copy()
        filterGene = self.db_.loc[self.db_.gene_name == geneExpr_].copy()
        givenDb = dbToMerge.merge(filterGene)
        XPop_ = patsy.dmatrix( self.formulaFE_ ,data = givenDb)
        YPop_ = givenDb.loc[:,self.expression_].copy()
        self.xMatrix = XPop_
        self.yVector = YPop_.values
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
    def calculateHerit(self,column_):
        totSum = 0
        lastRun = self.sigma2[len(self.sigma2)-1]
        for cov_ in list(lastRun):
            totSum += lastRun[cov_]
        self.h2 = lastRun[column_] / totSum
        self.finalDesiredSigma2 = lastRun[column_]
        self.finalTotalSigma = totSum
    def calculateAll(self,column_):
        finalTuple = []
        for geneExpr_ in self.genes_:
            self.defineFEM_PV(geneExpr_)
            self.defineREM(geneExpr_)
            if self.method_ == 'REML':
                self.estimateAddVarREML()
            elif self.method_ == 'ML':
                self.estimateAddVarML()
            else:
                print("Invalid option")
                break
            self.calculateHerit(column_)
            h2 = self.h2[0].copy()
            finalDesiredSigma2 = self.finalDesiredSigma2[0].copy()
            finalTotalSigma = self.finalTotalSigma[0].copy()
            pop = self.pop_
            maf = self.maf_
            hwe = self.hwe_
            vif = self.vif_
            sampInf = self.sampInf_
            namesRE = list(self.randomCovs)
            for num_,x in enumerate(namesRE):
                if num_ == 0:
                    RE = x
                else:
                    RE += "+" + x
            tuple_ = [geneExpr_,pop , maf , hwe , vif , sampInf,h2,finalDesiredSigma2,finalTotalSigma,self.method_,self.formulaFE_,RE]
            finalTuple.append(tuple_)
        finalDf = pd.DataFrame(finalTuple)
        finalDf.columns = ['Gene','Pop' , 'MAF' , 'HWE' , 'VIF' , 'sampleInference','HeritEst','DesiredVarEst','totalVarEst','method','formulaF','randEffects']
        self.results = finalDf
    def saveResults(self):
        check_ = 0
        try:
            df_ = pd.read_csv(self.saveRef_,sep = '|',header = 0)
        except:
            check_ = 1
        if check_ == 1:
            self.results.to_csv(self.saveRef_,index = False, header = True, sep = "|")
        else:
            df_ = df_.append(self.results)
            df_.drop_duplicates(inplace = True)
            df_.to_csv(self.saveRef_,index = False, header = True, sep = "|")



# Static params
pathPipelinePrograms = '/raid/genevol/users/lucas/heritability/01.Pipeline'
pathGCTA = '/raid/genevol/users/lucas/gcta'
pathTmp = '/scratch/genevol/users/lucas'
resultsFile = pathPipelinePrograms + '/Results/results.txt'
expression = 'tpm'
db = pd.read_csv('/raid/genevol/heritability/hla_expression.tsv',sep = '\t')
genes = sorted(db.gene_name.unique())
# Mutable params
samplesFile = pathPipelinePrograms + '/Samples/samples.filt' 
dbIndividuals = pd.read_csv(samplesFile,header = None,sep = '\t')
individuals = dbIndividuals.loc[:,[0]]
individuals.columns = ['subject_id']
pop = 'Geuvadis'
sampleInference = 'all'
maf = 0.01
hwe = 0.1
vif = 1.5
threads = 10
catCovs = ['sex','lab']
numCovs = None
extraRE = None
formulaFE = '~sex+lab'
parametersOpt = {'maxIt':50}

# Build Additive Variance Matrix - required to estimate narrow sense heritability
buildAVM = buildAdditiveVarianceMatrix(samplesFile=samplesFile,pop = pop,maf = maf , hwe = hwe, vif = vif,
                        threads=threads , pathPipeline=pathPipelinePrograms,
                        pathTmp=pathTmp , pathGCTA=pathGCTA
                        )
buildAVM.createTmpFolder()
# buildAVM.createBedFile()
buildAVM.createBedFileNPC()
buildAVM.calculateGCTA()
buildAVM.correctGRM()

# Narrow sense heritability calculation
# Calculating heritability for all required genes, in a classical manner
calculateH2GCTA = heritabilityGCTA(individuals = individuals,db = db , covs = catCovs , genes = genes,oldParams=buildAVM)
calculateH2GCTA.calculateAll()
# Paramaters needed to estimate in a more generalized way 
additiveMatrixList = dict()
buildAVM.readGRM()
additiveMatrixList['Genes'] = buildAVM.GRM
calculateH2Alt = heritabilityAlt(individuals = individuals,sampInf=sampleInference,db = db,formulaFE=formulaFE,expression = expression , genes = genes,oldParams=buildAVM,additiveMatrixList=additiveMatrixList,extraRE=extraRE,parametersOpt=parametersOpt,method='ML',resultsFile=resultsFile)
calculateH2Alt.calculateAll('Genes')
calculateH2Alt.saveResults()