import subprocess
import itertools
import os

def retStr(var_):
    if var_ != 'null':
        varStr_ = ''
        for element in var_:
            varStr_ += element + '_'
    else:
        varStr_ = var_
    return varStr_


# mafs = [0.01, 0.05]
mafs = [0.01]
# pops = [['FIN','GBR','TSI','CEU'],'null']
pops = [['FIN','GBR','TSI','CEU']]
vifs = [1.11]
# hwes = [1e-4,1e-7]
hwes = [1e-7]
sexs = ['null']
labs = ['null']
# outliers = [0.01 , 0.001 , 'null']
outliers = [0.01]
formulas = [{'fixed': '~pop + lab', 'random': 'null'}]




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
        f.write('formula: ' + str(formula_) + '\n')
    query_ = ['kedro', 'run' ,'--env=snpsParams_maf_' + maf_ + '_hwe_' + hwe_ + '_vif_' + vif_ + '_sampParams_' + popStr_ + 'sex_' + sexStr_ + '_lab_' + labStr_ + '_outliers_' + outliers_]
    print(query_)
    print(os.getcwd())
    subprocess.Popen(query_)

