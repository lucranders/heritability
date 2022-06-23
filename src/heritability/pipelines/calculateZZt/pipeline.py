"""
This is a boilerplate pipeline 'calculateZZt'
generated using Kedro 0.17.6
"""

from kedro.pipeline import pipeline, node,Pipeline
from .nodes import calculate_K_C_alpha, correct_K_C_alpha

def create_pipeline(**kwargs):
    templateK_C_alphaCalculus = Pipeline([
    node(
            func=calculate_K_C_alpha,
            inputs=['params:overrideSetsChrs',"selected_Sample","params:alphas","bed_Files"],
            outputs="Product_k_c_alpha",
            name="Calculate_product_k_c_alpha",
            tags=['Genotype','Matrix','Calculation']
    ),
    node(
            func=correct_K_C_alpha,
            inputs=['params:overrideSetsChrs',"selected_Sample","params:alphas","Product_k_c_alpha"],
            outputs="corrected_product_k_c_alpha",
            name="correction_product_k_c_alpha",
            tags=['Genotype','Matrix','Correction']
    ),
    ])
    setName_ = 'FullSet'
    fullSetChromosomes = pipeline(
        pipe=templateK_C_alphaCalculus,
        parameters={"params:overrideSetsChrs": "params:listChrsFull", "selected_Sample":"selected_Sample", "bed_Files":"bed_Files"},
        namespace="herit_" + setName_

    )
    setName_ = 'TwoSets'
    twoSetsChromosomes = pipeline(
        pipe=templateK_C_alphaCalculus,
        parameters={"params:overrideSetsChrs": "params:listChrsPartitioned2", "selected_Sample":"selected_Sample", "bed_Files":"bed_Files"},
        namespace="herit_" + setName_
    )

    setName_ = 'TwentyTwoSets'
    twentyTwoSetsChromosomes = pipeline(
        pipe=templateK_C_alphaCalculus,
        parameters={"params:overrideSetsChrs": "params:listChrsPartitioned22", "selected_Sample":"selected_Sample", "bed_Files":"bed_Files"},
        namespace="herit_" + setName_
    )

    entirepipeline = fullSetChromosomes + twoSetsChromosomes + twentyTwoSetsChromosomes
    return entirepipeline