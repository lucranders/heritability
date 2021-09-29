import subprocess
import os
import time
from datetime import datetime
from pathlib import Path
import pandas as pd

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

class pipelineGCTA:
    def __init__(self,pop,maf,hwe,vif,threads,genes,covs,db,pathPipeline,pathTmp,pathGCTA):
        self.pop_ = pop
        self.maf_ = maf
        self.hwe_ = hwe
        self.vif_ = vif
        self.genes_ = genes
        self.covs_ = covs
        self.threads_ = str(threads)
        self.pathPipeline_ = pathPipeline
        self.pathTmp_ = pathTmp
        self.pathGCTA_ = pathGCTA
        self.path_ = pathTmp + '/TempBed_pop_' + pop + "_maf_" + str(maf) + "_hwe_" + str(hwe) + "_vif_" + str(vif)
        self.db_ = db
    # Create .bed files
    def createTmpFolder(self):
        # Create temp folder
        subprocess.Popen(['mkdir',self.path_])
    def createFiles(self):
        basicCols = ['subject_id','subject_id']
        X_ = basicCols.copy()
        Y_ = basicCols.copy()
        Y_.append('tpm')
        print(self.db_.head(3))
        for cov_ in self.covs_:
            X_.append(cov_)
        if self.pop_ == 'NAfr':
            # Filtering known subjects information
            # Filtering only HLA-A because the same info repeats for each known gene expression
            XPop_ = self.db_.loc[(self.db_.pop != 'YRI')&(self.db_.gene_name == 'HLA-A'),X_]
            # Saving columns into temp folder
            XPop_.to_csv(self.path_ + '/X.txt',sep = ' ',index = False, header = False)
            for geneExpr_ in self.genes_:
                YPop_ = self.db_.loc[(self.db_.pop != 'YRI')&(self.db_.gene_name == geneExpr_),Y_]
                YPop_.to_csv(self.path_ + '/Y_' + geneExpr_.replace('-','') +'.txt',sep = ' ',index = False, header = False)
        else:
            filesamp = self.pathPipeline_ + '/Samples/samples.filt'
            filt = pd.read_csv(filesamp,sep = '\t',header=None)
            filt.columns = ['col1','col2']
            XPop_ = self.db_.loc[(self.db_.subject_id.isin(filt.col1))&(self.db_.gene_name == 'HLA-A'),X_]
            # Saving columns into temp folder
            XPop_.to_csv(self.path_ + '/X.txt',sep = ' ',index = False, header = False)
            for geneExpr_ in self.genes_:
                YPop_ = self.db_.loc[(self.db_.subject_id.isin(filt.col1))&(self.db_.gene_name == geneExpr_),Y_]
                YPop_.to_csv(self.path_ + '/Y_' + geneExpr_.replace('-','') +'.txt',sep = ' ',index = False, header = False)
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
        check_ = Path(pipeInit.path_ + '/statusGRMCorrection.txt')
        cond = check_.is_file()
        while not cond:
            cond = check_.is_file()
            time.sleep(1)
    def heritability(self):
        for geneExp_ in self.genes_:
            gctaLoc = '/raid/genevol/users/lucas/gcta/gcta64'
            method = '--reml'
            algotithm = '--reml-alg'
            algotithmV = '0'
            maxit = '--reml-maxit'
            maxitV = '1000'
            grm = '--grm'
            grmV = self.path_ + '/GCTA'
            pheno =  '--pheno'
            phenoV = self.path_ + '/Y_' + geneExp_.replace('-','') + '.txt'
            covs = '--covar'
            covsV = self.path_ + '/X.txt'
            fix = '--reml-est-fix'
            out = '--out'
            outV = self.path_ + '/Results' + geneExp_.replace('-','')
            cmd = [gctaLoc, method, algotithm, algotithmV, maxit, maxitV, grm, grmV, pheno, phenoV, covs, covsV, fix, out, outV]
            subprocess.Popen(cmd)



pathPipelinePrograms = '/raid/genevol/users/lucas/heritability/01.Pipeline'
pathGCTA = '/raid/genevol/users/lucas/gcta'
pathTmp = '/scratch/genevol/users/lucas'
pop = 'Geuvadis'
maf = 0.01
hwe = 1e-07
vif = 5
threads = 10
genes = ['HLA-A','HLA-B','HLA-C']
covs = ['sex','lab','pop']

db = pd.read_csv('/raid/genevol/heritability/hla_expression.tsv',sep = '\t')

pipeInit = pipelineGCTA(pop = pop,maf = maf , hwe = hwe, vif = vif,
                        genes = genes , covs = covs , db = db ,
                        threads=threads , pathPipeline=pathPipelinePrograms,
                        pathTmp=pathTmp , pathGCTA=pathGCTA
                        )

pipeInit.createTmpFolder()
pipeInit.createFiles()
pipeInit.createBedFile()
pipeInit.calculateGCTA()
pipeInit.correctGRM()
pipeInit.heritability()

# TODO:
    # parallelize process
