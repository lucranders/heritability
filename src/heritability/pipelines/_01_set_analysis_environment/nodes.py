from typing import Any, Dict
import subprocess
import os
import time
from datetime import datetime
import pandas as pd
import numpy as np
import pandas as pd
import logging
from kedro.config import ConfigLoader

file_name = 'logs/journals/set_analysis_environment.log'
logging.FileHandler(file_name)
logging.basicConfig(filename=file_name, encoding='utf-8', level=logging.DEBUG)

def retStr(var_):
    if var_ != None:
        varStr_ = ''
        for element in var_:
            varStr_ += element + '_'
    else:
        varStr_ = str(var_)
    return varStr_

logging.info('Current working directory:' + os. getcwd())
def createTempFolder(path_):
    # Create temporary folder
    conf_paths = ["conf/local"]
    conf_loader = ConfigLoader(conf_paths)
    parameters = conf_loader.get("paths*", "paths*/**")
    pathTemp = parameters['pathTemp']
    pathTempFiles = f'{pathTemp}/{path_}'
    proc_ = subprocess.Popen(['mkdir',pathTempFiles])
    proc_.wait()
    logging.info(f'Created temporary folder: {pathTempFiles}')
    return pathTempFiles

def select_sample(data: pd.DataFrame , genoData: pd.DataFrame, params_run: dict) -> Dict[str, Any]:
    """Node for selecting the desired sample to work with.
    The parameters are taken from conf/project/parameters.yml.
    The data and the parameters will be loaded and provided to this function
    automatically when the pipeline is executed and it is time to run this node.
    """
    pathTempFiles = createTempFolder(params_run['saveControl'])
    sampParams = params_run['sampParams']
    if sampParams['pop'] != None:
        data = data.loc[data['pop'].isin(sampParams['pop'])]
    if sampParams['sex'] != None:
        data = data.loc[data['sex'].isin(sampParams['sex'])]
    if sampParams['lab'] != None:
        data = data.loc[data['lab'].isin(sampParams['lab'])]
    # Merged data - intersection between genotyped and phenotyped individuals
    dfMerge = genoData.merge(data, how = 'inner')
    dfMerge.sort_values(['pop', 'subject_id'], inplace = True)
    sizeOriginalDf = dfMerge.drop_duplicates(['subject_id']).shape[0]
    logging.info('The original sample size is: ' + str(sizeOriginalDf))
    # Saving sample on temp folder
    if sampParams['outliers'] != None:
        if len(sampParams['genes']) > 1:
            from scipy.linalg import inv
            import numpy as np
            from scipy.stats import chi2
            cut_ = 1 - sampParams['outliers']
            gene_names = sampParams['genes']
            horizontalDf = dfMerge.pivot(index=['subject_id','sex','lab','pop'], columns=['gene_name'], values='tpm').reset_index()
            centeredValues = horizontalDf[gene_names] - horizontalDf[gene_names].mean()
            centeredValues.index = horizontalDf.subject_id
            covMatrix_ = np.cov(horizontalDf[gene_names].values.T)
            invCovMatrix_ = inv(covMatrix_)
            invCovMatrix_ = pd.DataFrame(invCovMatrix_)
            invCovMatrix_.index = gene_names
            invCovMatrix_.columns = gene_names
            mahalanobisMatrix = centeredValues @ invCovMatrix_ @ centeredValues.T
            mahalanobisFinal = np.diag(mahalanobisMatrix)
            mahalanobisFinal = pd.DataFrame(mahalanobisFinal)
            mahalanobisFinal.loc[:,'subject_id'] = horizontalDf.subject_id
            mahalanobisFinal.rename(columns = {0:'MahalanobisD'},inplace = True)
            mahalanobisFinal.loc[:,'Percentile'] = mahalanobisFinal.MahalanobisD.apply(lambda x: chi2.cdf(x,(len(gene_names))))
            mahalanobisFinal.loc[:,'cut'] = mahalanobisFinal.Percentile.apply(lambda x: 1 if x >= cut_ else 0)
            removed_sample_ = mahalanobisFinal.loc[mahalanobisFinal.cut == 1]
            logging.info(f'''
            Removed {removed_sample_.shape[0]} subjects:
            {removed_sample_}
            ''')
            dfSample = mahalanobisFinal.loc[mahalanobisFinal.cut == 0 ,['subject_id']].copy().drop_duplicates()
        else:
            df_gene_ = dfMerge.loc[dfMerge.gene_name == sampParams['genes'][0]]
            tab_ = df_gene_.describe()
            delta_ = tab_.loc['75%'] - tab_.loc['25%']
            upper_bound_ = tab_.loc['75%'] + 1.5*delta_
            lower_bound_ = tab_.loc['25%'] - 1.5*delta_
            chosen_sample_ = df_gene_.loc[(df_gene_.tpm >= lower_bound_.tpm)&(df_gene_.tpm <= upper_bound_.tpm)].copy()
            removed_sample_ = df_gene_.loc[~df_gene_.subject_id.isin(chosen_sample_.subject_id)]
            logging.info(f'''
            Gene: {sampParams['genes'][0]}
            Lower Boundary: {lower_bound_.tpm}
            Upper Boundary: {upper_bound_.tpm}
            Head: {chosen_sample_.head(3)}
            Removed {removed_sample_.shape[0]} subjects:
            {removed_sample_}
            ''')
            dfSample = chosen_sample_.loc[:, ['subject_id']].copy()
    else:
        dfSample = dfMerge.loc[:,['subject_id']].copy().drop_duplicates()
    sizeFinalDf = dfSample.shape[0]
    logging.info('The new sample size is: ' + str(sizeFinalDf))
    dfSample.loc[:,'1'] = dfSample.subject_id
    dfSample.to_csv(pathTempFiles + '/sample.txt', index = False, header= False, sep = ' ')
    final_ = {'selectedSample':dfSample.loc[:,['subject_id']] , 'sizeOriginalDf':sizeOriginalDf , 'sizeFinalDf':sizeFinalDf , 'originalDf': dfMerge, 'pathAnalysis': pathTempFiles}
    return final_

def filter_snps(params_run: dict, selectedSample: dict):
    pathTempFiles = selectedSample['pathAnalysis']
    snpsParams = params_run['snpsParams']
    conf_paths = ["conf/local"]
    conf_loader = ConfigLoader(conf_paths)
    parameters = conf_loader.get("paths*", "paths*/**")
    # Path with bad snps file -> indicates which snps are duplicate in all chromosomes
    pathDropSnps = parameters['pathDropSnps']
    pathPlink = parameters['plink']
    pathVcf = parameters['pathVcfFiles']
    main_query_ = pathPlink + " --vcf $input"
    if snpsParams['vif'] != None:
        main_query_ += " --indep 50 5 " + str(snpsParams['vif'])
    if snpsParams['maf'] != None:
        main_query_ += " --maf " + str(snpsParams['maf'])
    if snpsParams['hwe'] != None:
        main_query_ += " --hwe " + str(snpsParams['hwe'])
    if snpsParams['exclude'] != None:
        main_query_ += f" --exclude {pathDropSnps}"
    main_query_ += f" --keep {pathTempFiles}/sample.txt --mind 0.05 --geno 0.05 --vcf-half-call missing --out $filtered"
    # Create .bed files - 22 chromossomes
    query = f'''
        #!/bin/bash
        array_=()
        for chr_ in $(seq 1 22); do
            input={pathVcf}/ALL.chr"$chr_"_GRCh38.genotypes.20170504.vcf.gz
            filtered={pathTempFiles}/list_filt_snps_chr"$chr_"
            if test -f "$filtered.prune.in"; then
                echo "bed file for chromosome $chr_ already processed"
            else
                echo "creating bed file for chromosome $chr_"
                array_+=("{main_query_}")
            fi
        done;
        '''
    # Loop in batches of size 5 (parallel)
    query += '''
    let length_queries=${#array_[@]} step=5 how_many=length_queries/step
    let rest=length_queries-how_many*step

    for i in $(seq 0 $how_many); do
        let init=$i*step
        if (($i<$how_many)); then
            let end=($i+1)*$step-1
        else
            let end=$i*step+rest
        fi
        for batch_ in $(seq $init $end); do
            eval "${array_[$batch_]}"&
        done;
        wait;
    done;
        '''
    with open(f'{pathTempFiles}/prune_snps.sh', 'w') as f:
        f.write(query)
        f.close()
    cmd = ['chmod','+x', f'''{pathTempFiles}/prune_snps.sh''']
    p = subprocess.Popen(cmd)
    p.wait()
    cmd = ['bash',f'''{pathTempFiles}/prune_snps.sh''']
    logging.info(f'Initializing Snps pruning: bash script in {pathTempFiles}/prune_snps.sh')
    p = subprocess.Popen(cmd)
    p.wait()
    logging.info('Ended pruning step')
    return selectedSample

def remove_bad_snps(selectedSample):
    pathTempFiles = selectedSample['pathAnalysis']
    # Checking if all files were created
    filesPrunedSnps = [x for x in os.listdir(pathTempFiles) if 'list_filt_snps_chr' in x and '.prune.in' in x]
    logging.debug(f'''Pruned snps files: {filesPrunedSnps}''')
    for chr_ in range(1,23):
        prunedFile = f'list_filt_snps_chr{chr_}'
        snpsInfoFile = f'info_snps_{chr_}'
        # Reading pruned snps file
        dfPruned_ = pd.read_csv(f'{pathTempFiles}/{prunedFile}.prune.in', sep = '\t', header = None)
        dfPruned_.rename(columns = {0: 'SNP'}, inplace = True)
        # Reading snps info file
        dfSnpsInfo_ = pd.read_csv(f'data/01_raw/{snpsInfoFile}.txt', sep = '\t', header = None)
        logging.debug(f'''Comparing files {prunedFile} and {snpsInfoFile}''')
        dfSnpsInfo_.loc[:, 'aux_'] = 1
        count_duplicates = dfSnpsInfo_.groupby(['POS','SNP','CHR']).aux_.sum().reset_index()
        snps_ok_ = count_duplicates.loc[count_duplicates.aux_ == 1]
        snps_not_ok_ = count_duplicates.loc[count_duplicates.aux_ > 1]
        select_only_ = dfPruned_.loc[dfPruned_.SNP.isin(snps_ok_.SNP)].copy()
        logging.info(f'''
        There were {snps_not_ok_.shape[0]} snps with duplicates
        Main information:
            Chromosome: {chr_}
            Original snp count: {dfPruned_.shape[0]}
            Final snp count: {select_only_.shape[0]}
            Removed snps: {dfPruned_.shape[0] - select_only_.shape[0]}
        ''')
        select_only_.to_csv(f'{pathTempFiles}/{prunedFile}_dq.prune.in', header = None, index = False)
    return selectedSample


def create_bed_files(selectedSample):
    pathTempFiles = selectedSample['pathAnalysis']
    conf_paths = ["conf/local"]
    conf_loader = ConfigLoader(conf_paths)
    parameters = conf_loader.get("paths*", "paths*/**")
    pathPlink = parameters['plink']
    pathVcf = parameters['pathVcfFiles']
    mainQuery = f'''{pathPlink} --vcf $input --vcf-half-call missing --extract $filtered --keep {pathTempFiles}/sample.txt --make-bed --out $bed'''
    query = f'''#!/bin/bash
    array_=()
    for chr_ in $(seq 1 22); do
        input={pathVcf}/ALL.chr"$chr_"_GRCh38.genotypes.20170504.vcf.gz
        filtered={pathTempFiles}/list_filt_snps_chr"$chr_".prune.in
        bed={pathTempFiles}/pruned_chr$chr_
        if test -f "$bed.bed"; then
            echo "bed file for chromosome $chr_ already processed"
        else
            echo "creating bed file for chromosome $chr_"
            array_+=("{mainQuery}")
        fi
    done;

    '''
    # Loop in batches of size 5 (parallel)
    query += '''
    let length_queries=${#array_[@]} step=5 how_many=length_queries/step
    let rest=length_queries-how_many*step

    for i in $(seq 0 $how_many); do
        let init=$i*step
        if (($i<$how_many)); then
            let end=($i+1)*$step-1
        else
            let end=$i*step+rest
        fi
        for batch_ in $(seq $init $end); do
            eval "${array_[$batch_]}"&
        done;
        wait;
    done;
        '''
    with open(f'{pathTempFiles}/bedMaker.sh', 'w') as f:
        f.write(query)
        f.close()
    cmd = ['chmod','+x', f'''{pathTempFiles}/bedMaker.sh''']
    p = subprocess.Popen(cmd)
    p.wait()
    logging.info(f'Initializing creation of bed files from pruned snps: bash script in {pathTempFiles}/bedMaker.sh')
    cmd = ['bash',f'''{pathTempFiles}/bedMaker.sh''']
    p = subprocess.Popen(cmd)
    p.wait()
    return selectedSample