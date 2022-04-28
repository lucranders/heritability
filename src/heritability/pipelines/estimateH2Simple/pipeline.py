"""
This is a boilerplate pipeline 'estimateH2Simple'
generated using Kedro 0.17.6
"""

from kedro.pipeline import Pipeline, node
from .nodes import assembleParams , estimateSigmas2REMLSimple

def create_pipeline(**kwargs):
    return Pipeline([
        node(
            estimateSigmas2REMLSimple,
            ['params:snpsParams' , 'params:sampParams' , 'params:formula', 'params:GeneExpressions' , 'params:pathTempFiles',  'selected_Sample' ,'corrected_ZZt'],
            outputs='saved_sigma2_estimates',
            name="calculates_sigma2_given_different_initial_values",
        ),

    ])
