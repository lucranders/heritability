"""
This is a boilerplate pipeline 'estimateH2Simple'
generated using Kedro 0.17.6
"""

from kedro.pipeline import Pipeline, node
from .nodes import assembleParams , estimateSigmas2REMLSimple , estimateSigmas2MLSimple , estimateSigmas2MLSimpleR

def create_pipeline(**kwargs):
    return Pipeline([
        node(
            estimateSigmas2REMLSimple,
            ['params:snpsParams' , 'params:sampParams' , 'params:formula', 'params:GeneExpressions' , 'params:pathTempFiles',  'selected_Sample' ,'corrected_ZZt'],
            outputs='saved_sigma2_estimates_reml',
            name="calculates_sigma2_given_different_initial_values_reml",
        ),
        node(
            estimateSigmas2MLSimple,
            ['params:snpsParams' , 'params:sampParams' , 'params:formula', 'params:GeneExpressions' , 'params:pathTempFiles',  'selected_Sample' ,'corrected_ZZt'],
            outputs='saved_sigma2_estimates_ml',
            name="calculates_sigma2_given_different_initial_values_ml",
        ),
        # node(
        #     estimateSigmas2MLSimpleR,
        #     ['params:snpsParams' , 'params:sampParams' , 'params:formula', 'params:GeneExpressions' , 'params:pathTempFiles',  'selected_Sample' ,'corrected_ZZt'],
        #     outputs='saved_sigma2_estimates_ml_r',
        #     name="calculates_sigma2_given_different_initial_values_ml_r",
        # ),
    ])