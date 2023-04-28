import subprocess
import itertools
import os
from kedro.config import ConfigLoader
import time
from os.path import exists
import yaml



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
            

mafs = [0.05,.1]
# mafs = [0.01]
pops = [['FIN','GBR','TSI','CEU','YRI'], ['FIN','GBR','TSI','CEU']]
# pops = [['FIN','GBR','TSI','CEU','YRI']]
vifs = [1.5,2]
hwes = [1e-7,1e-4]
# hwes = [1e-7]
sexs = [None]
labs = [None]
genes = [
        ["HLA-A", "HLA-B", "HLA-C", "HLA-DQA1", "HLA-DPB1", "HLA-DRB1", "HLA-DPA1", "HLA-DQB1", "HLA-DRA"],
        ['HLA-A']
        ]

outliers = [0.01, None]
formulas = [{'fixed': ['pop+lab'], 'random': None}]
matrices_ = {
                # 'K_C':[-1],
                # 'std_K_C':[-1]
                'std_GCTA': 10,
                'GCTA': 10
            }


sets_ = {
        'pipelineFullChrSer': {'all': [x for x in range(1,23)]}
        }



for element in itertools.product(mafs , pops , vifs , hwes , sexs , labs , outliers,genes,formulas):
    maf_ = element[0]
    pop_ = element[1]
    vif_ = element[2]
    hwe_ = element[3]
    sex_ = element[4]
    lab_ = element[5]
    outliers_ = element[6]
    genes_ = element[7]
    formula_ = element[8]
    print(element)
    popStr_ = "_".join(pop_)
    sexStr_ = "_".join(sex_) if sex_ is not None else 'null'
    labStr_ = "_".join(lab_) if lab_ is not None else 'null'
    geneStr_ = "_".join(genes_)
    folder_name_ = f'''snpsParams_maf_{maf_}_hwe_{hwe_}_vif_{vif_}_sampParams_{popStr_}_sex_{sexStr_}_lab_{labStr_}_outliers_{'null' if outliers_ is None else outliers_}_genes_{geneStr_}'''
    path_ =  f'''conf/{folder_name_}''' 
    if not exists(path_):
        proc_ = subprocess.Popen(['mkdir',path_])
        proc_.wait()
    else:
        print('env already created!\n' + path_)
    final_dict_ = {
                    'params_run':{
                                'snpsParams': {
                                    'maf': maf_,
                                    'hwe': hwe_,
                                    'vif': vif_,
                                    'exclude': 'yes',
                                    },
                                'sampParams':{
                                    'pop': pop_,
                                    'lab': lab_,
                                    'sex': sex_,
                                    'outliers': outliers_,
                                    'genes': genes_
                                },
                                'formula': formula_,
                                'saveControl': folder_name_,
                                'matrices': {
                                            'sets_': sets_, 
                                            'type_': matrices_
                                            }
                                }
                }
    with open(path_+'/parameters.yml','w') as file:
        yaml.dump(final_dict_, file,sort_keys=True)
        file.close()
    query_ = ['kedro', 'run' ,f'--env={folder_name_}']
    # '--to-nodes=_4:matrices_calculation.calculate_matrices']
    print(query_)
    proc_ = subprocess.Popen(query_)
    proc_.wait()



