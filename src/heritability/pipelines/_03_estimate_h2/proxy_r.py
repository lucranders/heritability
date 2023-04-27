import subprocess
import logging

logging.basicConfig(filename='logs/journals/estimate_h2.log', encoding='utf-8', level=logging.DEBUG)


def estimate_herit_R(params, selectedSample):
    pathAnalysis = selectedSample['pathAnalysis']
    logging.info("Initializing heritability estimates step")
    matrices_ = []
    for matrix_type_ in params['matrices']['type_']:
        if type(params['matrices']['type_'][matrix_type_]) == list:
            for element_ in params['matrices']['type_'][matrix_type_]:
                print(element_)
                matrices_.append(f'''{matrix_type_}_{element_}''')
        else:
            matrices_.append(f'''{matrix_type_}''')
    query = f'''
                #!/bin/bash
                array_=()
                path_infos={pathAnalysis}
                formulas=({' '.join(params['formula']['fixed'])})
                sets=({' '.join(params['matrices']['sets_'])})
                matrices=({' '.join(matrices_)})
                genes=({' '.join(params['sampParams']['genes'])})
    '''
    query+='''
                for formula_ in ${formulas[@]}; do
                    for set_ in ${sets[@]}; do
                        for matrix_ in ${matrices[@]}; do
                            for gene_ in ${genes[@]}; do
                                name_file=$path_infos/resultfit_$formula_$set_$matrix_$gene.txt
                                if test -f "$name_file"; then
                                    echo "$name_file already exists and will not be processed"
                                else
                                    echo "$name_file still not exists and will be processed"
                                    array_+=("Rscript src/heritability/pipelines/_03_estimate_h2/calculate_herit.r $path_infos $formula_ $set_ $matrix_ $gene_")
                                fi
                            done;
                        done;
                    done;
                done;

                let length_queries=${#array_[@]} step=10 how_many=length_queries/step
                let rest=length_queries-how_many*step
                for i in $(seq 0 $how_many); do
                    let init=$i*step
                    if (($i<$how_many)); then
                        let end=($i+1)*$step-1
                    else
                        let end=$i*step+rest
                    fi
                    for batch_ in $(seq $init $end); do
                        echo ${array_[$batch_]}
                        eval "${array_[$batch_]}"&
                    done;
                    wait;
                done;
    '''
    with open(f'{pathAnalysis}/estimate_heritability.sh', 'w') as f:
        f.write(query)
        f.close()
    cmd = ['chmod','+x', f'''{pathAnalysis}/estimate_heritability.sh''']
    p = subprocess.Popen(cmd)
    p.wait()
    cmd = ['bash',f'''{pathAnalysis}/estimate_heritability.sh''']
    logging.info(f'Initializing Snps pruning: bash script in {pathAnalysis}/estimate_heritability.sh')
    p = subprocess.Popen(cmd)
    p.wait()
    logging.info('Ended pruning step')
    return 1