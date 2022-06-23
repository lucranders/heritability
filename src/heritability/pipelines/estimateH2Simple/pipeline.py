"""
This is a boilerplate pipeline 'estimateH2Simple'
generated using Kedro 0.17.6
"""

from kedro.pipeline import pipeline, node,Pipeline
from .nodes import assembleParams , estimateSigmas2REMLSimple, estimateSigmas2REMLSimpleSingleStart,execAllSimple

def create_pipeline(**kwargs):
    templateHeritReml = Pipeline([
    node(
            func=estimateSigmas2REMLSimpleSingleStart,
            inputs=['params:snpsParams' , 'params:sampParams' , 'params:formula', 'params:GeneExpressions',  'selected_Sample' ,'params:alphas','corrected_product_k_c_alpha','params:overrideSetsChrs'],
            outputs="heritability_estimates_for_different_alphas",
            name="calculate_heritability_estimates_for_different_alphas",
            tags=['Heritability','Estimation']
    ),
    ])
    setName1_ = 'FullSet'
    fullSetChromosomes = pipeline(
        pipe=templateHeritReml,
        inputs={'selected_Sample':'selected_Sample','corrected_product_k_c_alpha':"herit_" + setName1_+'.corrected_product_k_c_alpha'},
        parameters={"params:overrideSetsChrs": "params:listChrsFull"},
        namespace="herit_" + setName1_
    )
    setName2_ = 'TwoSets'
    twoSetsChromosomes = pipeline(
        pipe=templateHeritReml,
        inputs={'selected_Sample':'selected_Sample','corrected_product_k_c_alpha':"herit_" + setName2_+'.corrected_product_k_c_alpha'},
        parameters={"params:overrideSetsChrs": "params:listChrsFull"},
        namespace="herit_" + setName2_
    )
    setName3_ = 'TwentyTwoSets'
    twentyTwoSetsChromosomes = pipeline(
        pipe=templateHeritReml,
        inputs={'selected_Sample':'selected_Sample','corrected_product_k_c_alpha':"herit_" + setName3_+'.corrected_product_k_c_alpha'},
        parameters={"params:overrideSetsChrs": "params:listChrsFull"},
        namespace="herit_" + setName3_
    )
    endPipe = Pipeline([
    node(
            func=execAllSimple,
            inputs=["herit_" + setName1_ + '.heritability_estimates_for_different_alphas',"herit_" + setName2_ + '.heritability_estimates_for_different_alphas',"herit_" + setName3_ + '.heritability_estimates_for_different_alphas','params:saveControl'],
            outputs="saved_Files",
            name="saves_all_files",
            tags=['Heritability','Estimation'],
            namespace='save_file'
    ),
    ])
    entirepipeline = fullSetChromosomes + twoSetsChromosomes + twentyTwoSetsChromosomes + endPipe
    return entirepipeline




