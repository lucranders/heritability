"""
This is a boilerplate pipeline 'estimateH2Simple'
generated using Kedro 0.17.6
"""

from kedro.pipeline import pipeline, node,Pipeline
from .proxy_r import estimate_herit_R
# from .nodes import assembleParams , estimateSigmas2REMLSimple, estimateSigmas2REMLSimpleSingleStart,estimateSigmas2REMLSimpleSingleStartGCTA

def create_pipeline(**kwargs):
    templateHeritReml = Pipeline([
    node(
            func=estimate_herit_R,
            inputs=['params:params_run',"input_data"],
            outputs="heritability_estimates",
            name="calculate_heritability_estimates",
            tags=['Heritability','Estimation']
    )
    ])
    estimate_herit_pipe = pipeline(
	pipe=templateHeritReml,
	inputs={'input_data':'_4:matrices_calculation.non_negative_specified_matrices'},
	namespace="_5:heritability_calculation"
	)
    return estimate_herit_pipe




