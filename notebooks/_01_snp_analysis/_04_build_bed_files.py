import subprocess
import os
from kedro.config import ConfigLoader

# Setting paths
conf_paths = ["../../conf/local"]
conf_loader = ConfigLoader(conf_paths)
parameters = conf_loader.get("paths*", "paths*/**")
pathPlink = parameters['plink']
pathVcf = parameters['pathVcfFiles']
pathTemp = parameters['pathTemp']
for name_ in ['test_0_0','test_0_1']:
    pathTempFiles = f'''{pathTemp}/{name_}'''
    # Query to produce new bed files ->
    # From original .prune.in files and from
    # the new list (drop duplicates) .prune.in file
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
            let end=($i+1)*$step
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
    # cmd = ['chmod','+x', f'''{pathTempFiles}/bedMaker.sh''']
    # p = subprocess.Popen(cmd)
    # p.wait()
    cmd = ['bash',f'''{pathTempFiles}/bedMaker.sh''']
    p = subprocess.Popen(cmd)
    p.wait()


