import subprocess
import itertools
import os
from kedro.config import ConfigLoader
import time

def retStr(var_):
    if var_ != 'null':
        varStr_ = ''
        for element in var_:
            varStr_ += element + '_'
    else:
        varStr_ = var_
    return varStr_
def checkLogSizes(path_,name_,ext_):
    listSizes = []
    for chr_ in range(1,23):
        try:
            size = os.path.getsize(path_ + '/' + name_ + str(chr_) + ext_)
        except:
            size = 0
    
        listSizes.append(size > 0)
    return listSizes
def wait(pathTempFiles,name_,ext_):
    sizes = checkLogSizes(pathTempFiles,name_,ext_)
    cond = sum(sizes) < 22
    # If not, waits until so
    while cond:
        sizes = checkLogSizes(pathTempFiles,name_,ext_)
        cond = sum(sizes) < 22
        if cond:
            # keep waiting to launch next request
            time.sleep(10)
            print("Waiting")
        else:
            print("Done!")
            
mafs = [0.01, 0.05]
# mafs = [0.01]
pops = [['FIN','GBR','TSI','CEU'],'null']
# pops = [['FIN','GBR','TSI','CEU']]
vifs = [1.11]
hwes = [1e-4,1e-7]
# hwes = [1e-7]
sexs = ['null']
labs = ['null']
outliers = [0.01 , 0.001 , 'null']
formulas = [{'fixed': '~pop + lab', 'random': None}]




for element in itertools.product(mafs , pops , vifs , hwes , sexs , labs , outliers,formulas):
    maf_ = str(element[0])
    pop_ = element[1]
    vif_ = str(element[2])
    hwe_ = str(element[3])
    sex_ = element[4]
    lab_ = element[5]
    outliers_ = str(element[6])
    formula_ = element[7]
    print(element)
    popStr_ = retStr(pop_)
    sexStr_ = retStr(sex_)
    labStr_ = retStr(lab_)
    path_ =  'conf/snpsParams_maf_' + maf_ + '_hwe_' + hwe_ + '_vif_' + vif_ + '_sampParams_' + popStr_ + 'sex_' + sexStr_ + '_lab_' + labStr_ + '_outliers_' + outliers_
    proc_ = subprocess.Popen(['mkdir',path_])
    proc_.wait()
    with open(path_+'/parameters.yml','w') as f:
        f.write('# Parameters for snp selection\n')
        f.write("snpsParams: {'maf': " + maf_ + ", 'hwe': " + hwe_ + ",'vif': " + vif_ + '}\n')
        f.write('# Sample parameters\n')
        f.write("sampParams: {'pop': " + str(pop_) + ", 'sex': " + str(sex_) + ",'lab': " + str(lab_) + ", 'outliers': " + outliers_ + '}\n')
        f.write('# Chosen model\n')
        f.write('formula: ' + str(formula_).replace('None','Null') + '\n')
    query_ = ['kedro', 'run' ,'--env=snpsParams_maf_' + maf_ + '_hwe_' + hwe_ + '_vif_' + vif_ + '_sampParams_' + popStr_ + 'sex_' + sexStr_ + '_lab_' + labStr_ + '_outliers_' + outliers_]
    print(query_)
    print(os.getcwd())
    subprocess.Popen(query_)
    conf_paths = ["conf/local"]
    conf_loader = ConfigLoader(conf_paths)
    parameters = conf_loader.get("paths*", "paths*/**")
    pathTemp = parameters['pathTemp']
    pathTempFiles = pathTemp + '/Temp_snpsParams_maf_' + maf_ + '_hwe_' + hwe_ + '_vif_' + vif_ + '_sampParams_' + popStr_ + 'sex_' + sexStr_ + '_lab_' + labStr_ + '_outliers_' + outliers_
    wait(pathTempFiles,'chr','.bed')