import subprocess

pop = 'NAfr'
maf = 0.01
hwe = 1e-08
vif = 5
cmd = ['qsub', '-v' ,'pop=' + pop + ',maf_=' + str(maf) +',hwe_=' + str(hwe) + ',vif_=' + str(vif) ,'/raid/genevol/users/lucas/heritability/01.Pipeline/01.DataPrepGeneral.sh']
subprocess.Popen(cmd).wait()

# TODO: Create temporary folder and then keep checking if '.log' archives are all of size > 0 (Or lack of 'temporary bed' files). If so, proceed to GCTA