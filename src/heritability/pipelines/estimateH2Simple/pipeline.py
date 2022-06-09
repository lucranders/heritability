"""
This is a boilerplate pipeline 'estimateH2Simple'
generated using Kedro 0.17.6
"""

from kedro.pipeline import Pipeline, node
from .nodes import assembleParams , estimateSigmas2REMLSimple, estimateSigmas2REMLSimpleSingleStart,execAllSimple

def create_pipeline(**kwargs):
    return Pipeline([
        # node(
        #     estimateSigmas2REMLSimple,
        #     ['params:snpsParams' , 'params:sampParams' , 'params:formula', 'params:GeneExpressions',  'selected_Sample' ,'corrected_ZZt' , 'params:saveControl'],
        #     outputs='saved_sigma2_estimates_reml',
        #     name="calculates_sigma2_given_different_initial_values_reml",
        # ),
        node(
            estimateSigmas2REMLSimpleSingleStart,
            ['params:snpsParams' , 'params:sampParams' , 'params:formula', 'params:GeneExpressions',  'selected_Sample' ,'corrected_ZZt'],
            outputs='saved_sigma2_estimates_reml',
            name="calculates_1_gene_components_reml",
        ),
        node(
            estimateSigmas2REMLSimpleSingleStart,
            ['params:snpsParams' , 'params:sampParams' , 'params:formula', 'params:GeneExpressions',  'selected_Sample' ,'corrected_ZZt_partitioned_2_components'],
            outputs='saved_sigma2_estimates_reml_2_components',
            name="calculates_2_gene_components_reml",
        ),
        node(
            estimateSigmas2REMLSimpleSingleStart,
            ['params:snpsParams' , 'params:sampParams' , 'params:formula', 'params:GeneExpressions',  'selected_Sample' ,'corrected_ZZt_partitioned_22_components'],
            outputs='saved_sigma2_estimates_reml_22_components',
            name="calculates_22_gene_components_reml",
        ),
        node(
            execAllSimple,
            ['saved_sigma2_estimates_reml','saved_sigma2_estimates_reml_2_components','saved_sigma2_estimates_reml_22_components','params:saveControl'],
            outputs='all_estimates_from_additive_models',
            name="calculates_all_estimates_from_additive_models",
        ),
    ])

