"""
This is a boilerplate pipeline 'estimateH2Simple'
generated using Kedro 0.17.6
"""
import .estimateH2

def create_class_herit_calc(snpsParams: dict , sampParams: dict , formula: dict):
    vif_ = snpsParams['vif']
    maf_ = snpsParams['maf']
    hwe_ = snpsParams['hwe']
    pop_ = ''
    for i,x in enumerate( sampParams['pop'] ):
        if i < ( len(sampParams['pop']) - 1 ):
            pop_ += x + '_'
        else:
            pop_ += x
    outLiers = sampParams['outliers']
        
def define_FEM_PV()