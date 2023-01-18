import subprocess
import os
from kedro.config import ConfigLoader

conf_paths = ["../../conf/local"]
conf_loader = ConfigLoader(conf_paths)
parameters = conf_loader.get("paths*", "paths*/**")
name_ = 'test_'
pathPlink = parameters['plink']
pathVcf = parameters['pathVcfFiles']
pathTemp = parameters['pathTemp']
pathTempFiles = f'''{pathTemp}/{name_}'''
mainQuery = f'''{pathPlink} --vcf $input --vcf-half-call missing --extract $filtered --keep {pathTempFiles}/sample.txt --make-bed --out $bed'''
filesTemp = [x for x in os.listdir(pathTempFiles)]
for chr_ in range(1,23):
    stdPruned = f'list_filt_snps_chr{chr_}_test_.prune.in'
    dedupPruned = f'new_list_filt_snps_{chr_}.prune.in'
    if f'stdPruned_chr{chr_}.bed' not in filesTemp:
        stdQuery = f'''
input={pathVcf}/ALL.chr{chr_}_GRCh38.genotypes.20170504.vcf.gz
filtered={stdPruned}
bed={pathTempFiles}/stdPruned_chr{chr_}
eval {mainQuery}
                        '''
        print(stdQuery)
        with open(pathTempFiles + '/stdPrunedBedMaker' + str(chr_) + '.sh', 'w') as f:
            f.write(stdQuery)
            f.write('\n')
            f.close()
        cmd = ['sh', pathTempFiles + '/stdPrunedBedMaker' + str(chr_) + '.sh']
        subprocess.Popen(cmd)
    if f'dedupPruned_chr{chr_}.bed' not in filesTemp:
        dedupQuery = f'''
input={pathVcf}/ALL.chr{chr_}_GRCh38.genotypes.20170504.vcf.gz
filtered={dedupPruned}
bed={pathTempFiles}/dedupPruned_chr{chr_}
eval {mainQuery}
                        '''
        print(dedupQuery)
        with open(pathTempFiles + '/dedupPrunedBedMaker' + str(chr_) + '.sh', 'w') as f:
            f.write(dedupQuery)
            f.write('\n')
            f.close()
        cmd = ['sh', pathTempFiles + '/dedupPrunedBedMaker' + str(chr_) + '.sh']
        p = subprocess.Popen(cmd)
        p.wait()

