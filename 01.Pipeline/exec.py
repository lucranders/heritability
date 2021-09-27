import subprocess
import os
import time
from datetime import datetime

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
def updateLog(path_,status_):
    
    if status_ == 0:
        msg_ = str(datetime.now()) + ": Waiting..."
    else:
        msg_ = str(datetime.now()) + ": Done!"
    
    f = open(path_ + "/status.txt", "a")
    f.write(msg_ + '\n')
    f.close()

# Create .bed files
def createBedFile(pop_, maf_ , hwe_ , vif_):

    cmd = ['qsub', '-v' ,'pop=' + pop_ + ',maf_=' + str(maf_) +',hwe_=' + str(hwe_) + ',vif_=' + str(vif_) ,'/raid/genevol/users/lucas/heritability/01.Pipeline/01.DataPrepGeneral.sh']
    subprocess.Popen(cmd)

def calculateGCTA(path_,pop_,threads_):
    
    if pop_ == 'NAfr':
        filesamp='/raid/genevol/users/lucas/heritability/01.Pipeline/Samples/nAfr.filt'
        cmd = ['/raid/genevol/users/lucas/gcta/gcta64', '--mbfile' ,path_ + '/chrs.txt','--keep',filesamp,'--make-grm','--out','GCTA','--thread-num',threads_]
    elif pop_ == 'Geuvadis':
        filesamp='/raid/genevol/users/lucas/heritability/01.Pipeline/Samples/samples.filt'
        cmd = ['/raid/genevol/users/lucas/gcta/gcta64', '--mbfile' ,path_ + '/chrs.txt','--keep',filesamp,'--make-grm','--out','GCTA','--thread-num',threads_]
    else:
        cmd = ['/raid/genevol/users/lucas/gcta/gcta64', '--mbfile' ,path_ + '/chrs.txt','--make-grm','--out','GCTA','--thread-num',threads_]
    
    subprocess.Popen(cmd)


# Create GRM, after creation of .bed file
def createGRM(pop_, maf_ , hwe_ , vif_ , threads_):
    
    # Create temp folder and .bed files
    createBedFile(pop_, maf_ , hwe_ , vif_)
    path_ = '/scratch/genevol/users/lucas/TempBed_pop_' + pop_ + "_maf_" + str(maf_) + "_hwe_" + str(hwe_) + "_vif_" + str(vif_)

    # Compute amount of chromossomes processed
    sizes = checkLogSizes(path_)
    cond = sum(sizes) < 22

    # Wait until createBedFile is processed
    while cond:
        sizes = checkLogSizes(path_)
        cond = sum(sizes) < 22
        if cond:
            updateLog(path_,0)
            time.sleep(10)
        else:
            updateLog(path_,1)

    for chr_ in range(1,23):
        f = open(path_ + '/chrs.txt', "a")
        f.write('chr' + str(chr_) + '\n')
        f.close()
    
    calculateGCTA(path_,pop_,str(threads_))




pop = 'NAfr'
maf = 0.01
hwe = 1e-08
vif = 5
threads = 10

createGRM(pop,maf,hwe,vif,threads)

# TODO:
    # Waiting process to GCTA calculation
    # Calculate heritability
    # parallelize process