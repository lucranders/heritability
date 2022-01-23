import subprocess
import os
import time
from datetime import datetime
from pathlib import Path
import pandas as pd
import numpy as np
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri


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
    def __init__(self,pop,maf,hwe,vif,threads,pathPipeline,pathTmp,pathGCTA):
        self.pop_ = pop
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
    def createBedFile(self):
        # Create .bed files - 22 chromossomes
        cmd = ['qsub', '-v' ,'pop=' + self.pop_ + ',maf_=' + str(self.maf_) +',hwe_=' + str(self.hwe_) + ',vif_=' + str(self.vif_) , self.pathPipeline_ + '/01.DataPrepGeneral.sh']
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
    def calculateGCTA(self):
        if self.pop_ == 'NAfr':
            filesamp = self.pathPipeline_ + '/Samples/nAfr.filt'
            cmd = [self.pathGCTA_ + '/gcta64', '--mbfile' ,self.path_ + '/chrs.txt','--keep',filesamp,'--make-grm','--out',self.path_+'/GCTA','--thread-num',self.threads_]
        elif self.pop_ == 'Geuvadis':
            filesamp = self.pathPipeline_ + '/Samples/samples.filt'
            cmd = [self.pathGCTA_ + '/gcta64', '--mbfile' ,self.path_ + '/chrs.txt','--keep',filesamp,'--make-grm','--out',self.path_+'/GCTA','--thread-num',self.threads_]
        else:
            cmd = [self.pathGCTA_ + '/gcta64', '--mbfile' ,self.path_ + '/chrs.txt','--make-grm','--out',self.path_+'/GCTA','--thread-num',self.threads_]
        subprocess.Popen(cmd)
        # check whether GRM binaries already exists - means process is finished
        check_ = Path(self.path_ + '/GCTA.grm.bin')
        cond = check_.is_file()
        while not cond:
            cond = check_.is_file()
            if not cond:
                # Updates status (inside tmp folder)
                updateLog(self.path_,0,'GRMStatus.txt')
                time.sleep(10)
            else:
                # Updates status (inside tmp folder)
                updateLog(self.path_,1,'GRMStatus.txt')
    def correctGRM(self):
        subprocess.Popen(['Rscript',self.pathPipeline_ + '/02.matrixCorrection.r',self.path_])
        check_ = Path(self.path_ + '/statusGRMCorrection.txt')
        cond = check_.is_file()
        while not cond:
            cond = check_.is_file()
            time.sleep(1)
        subprocess.Popen(['rm',self.path_ + '/statusGRMCorrection.txt'])
        pandas2ri.activate()
        readRDS = robjects.r['readRDS']
        self.correctedGRM = readRDS(self.path_ + '/GCTA.rds')

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
    def defineREM(self,currentGene):
        self.randomCovs = dict()
        self.randomCovs['Genes'] = self.correctedGRM.copy()
        if self.extraRE != None:
            for var_ in self.extraRE:
                oneHotEncoding = pd.get_dummies(self.db_.loc[self.db_.gene_name == currentGene , var_]).values
                covMat = np.matmul( oneHotEncoding , np.transpose(oneHotEncoding) )
                self.randomCovs[var_] = covMat.copy()
                del(covMat)
        self.randomCovs['Residuals'] = np.diag(np.full(self.correctedGRM.shape[0],1))






pathPipelinePrograms = '/raid/genevol/users/lucas/heritability/01.Pipeline'
pathGCTA = '/raid/genevol/users/lucas/gcta'
pathTmp = '/scratch/genevol/users/lucas'
samplesFile = pathPipelinePrograms + '/Samples/samples.filt' 
pop = 'Geuvadis'
maf = 0.01
hwe = 1e-07
vif = 5
threads = 10
genes = ['HLA-A','HLA-B','HLA-C']
covs = ['sex','lab','pop']
dbIndividuals = pd.read_csv(samplesFile,header = None,sep = '\t')
individuals = dbIndividuals.loc[:,[0]]
individuals.columns = ['subject_id']
db = pd.read_csv('/raid/genevol/heritability/hla_expression.tsv',sep = '\t')

buildAVM = buildAdditiveVarianceMatrix(pop = pop,maf = maf , hwe = hwe, vif = vif,
                        threads=threads , pathPipeline=pathPipelinePrograms,
                        pathTmp=pathTmp , pathGCTA=pathGCTA
                        )
calculateH2GCTA = heritabilityGCTA(individuals = individuals,db = db , covs = covs , genes = genes,oldParams=buildAVM)
buildAVM.createTmpFolder()
buildAVM.createBedFile()
buildAVM.calculateGCTA()
buildAVM.correctGRM()

calculateH2GCTA.calculateAll()

# TODO:
    # parallelize process
