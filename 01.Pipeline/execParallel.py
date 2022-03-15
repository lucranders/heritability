from distutils.command.build import build
import heritabilityFunctions as herit
import pandas as pd
from multiprocessing import Pool, freeze_support

# Paths
pathPipelinePrograms = '/raid/genevol/users/lucas/heritability/01.Pipeline'
pathGCTA = '/raid/genevol/users/lucas/gcta'
pathTmp = '/scratch/genevol/users/lucas'
resultsFile = pathPipelinePrograms + '/Results/resultsNew.txt'

# Static Params
# Variable name (transcriptions per million)
expression = 'tpm'
# dataframe, containing gene expression information + subjects information
db = pd.read_csv('/raid/genevol/heritability/hla_expression.tsv',sep = '\t')
# Genes names (whose expressions were measured) 
genes = sorted(db.gene_name.unique())
# Individuals to filter
individualsDf = pd.read_csv(pathPipelinePrograms + '/Samples/europeans.txt', header = None, sep = ' ')
individualsDf.columns = ['subject_id']
individuals = list(individualsDf.loc[:,'subject_id'])


# Mutable params
pop = 'European'
sampleInference = 'European'
mafs = [0.01, 0.05]
hwes = [1e-7, 1e-4]
vifs = [1.11]
args_ = []
for maf_ in mafs:
    for hwe_ in hwes:
        for vif_ in vifs:
            args_.append((maf_,hwe_,vif_))
threads = 10
catCovs = None
numCovs = None
extraRE = None
formulaFE = '~sex+lab'
parametersOpt = {'maxIt':100,'conv':1e-4}

def parallel(maf,hwe,vif):
    # Build Additive Variance Matrix - required to estimate narrow sense heritability
    buildAVM = herit.buildAdditiveVarianceMatrix(sample=individuals,pop = pop,maf = maf , hwe = hwe, vif = vif,
                            threads=threads , pathPipeline=pathPipelinePrograms,
                            pathTmp=pathTmp , pathGCTA=pathGCTA
                            )
    buildAVM.createTmpFolder()
    # buildAVM.createBedFile()
    buildAVM.createBedFileNPC()
    buildAVM.createChrRef()
    buildAVM.calculateGCTA()
    buildAVM.correctGRM()
    # Narrow sense heritability calculation
    # Calculating heritability for all required genes, in a classical manner
    calculateH2GCTA = herit.heritabilityGCTA(individuals = individualsDf,db = db , covs = ['sex','lab'] , genes = genes,oldParams=buildAVM)
    calculateH2GCTA.calculateAll(suffix_='_correction')
    # Paramaters needed to estimate in a more generalized way 
    additiveMatrixList = dict()
    buildAVM.readGRM()
    additiveMatrixList['Genes'] = buildAVM.GRM
    method = 'ML'
    calculateH2Alt = herit.heritabilityAlt(individuals = individualsDf,sampInf=sampleInference,db = db,formulaFE=formulaFE,expression = expression , genes = genes,oldParams=buildAVM,additiveMatrixList=additiveMatrixList,extraRE=extraRE,parametersOpt=parametersOpt,method=method,resultsFile=resultsFile)
    calculateH2Alt.calculateAll('Genes',compareGCTA=True)
    calculateH2Alt.saveResults()
    method = 'REML'
    calculateH2AltREML = herit.heritabilityAlt(individuals = individualsDf,sampInf=sampleInference,db = db,formulaFE=formulaFE,expression = expression , genes = genes,oldParams=buildAVM,additiveMatrixList=additiveMatrixList,extraRE=extraRE,parametersOpt=parametersOpt,method=method,resultsFile=resultsFile)
    calculateH2AltREML.calculateAll('Genes',compareGCTA=True)
    calculateH2AltREML.saveResults()

freeze_support()
with Pool() as pool:
    run = pool.starmap(parallel, args_)