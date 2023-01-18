"""Project pipelines."""
from typing import Dict

from kedro.pipeline import Pipeline

from heritability.pipelines import setAnalysisEnvironment as sae
from heritability.pipelines import calculateZZt as cZZt
from heritability.pipelines import estimateH2Simple as eH2S


def register_pipelines() -> Dict[str, Pipeline]:
    """Register the project's pipelines.

    Returns:
        A mapping from a pipeline name to a ``Pipeline`` object.

    """
    setAnalysisEnvironmentPipeline = sae.create_pipeline()
    calculateZZt = cZZt.create_pipeline()
    estimateH2Simple = eH2S.create_pipeline()

    setName_ = 'FullSet'
    pipelineFullChrSer = setAnalysisEnvironmentPipeline + \
                            calculateZZt.only_nodes_with_namespace("herit_" + setName_) +\
                            estimateH2Simple.only_nodes_with_namespace("herit_" + setName_)
    setName_ = 'TwoSets'
    twoSetsChromosomes = setAnalysisEnvironmentPipeline + \
                            calculateZZt.only_nodes_with_namespace("herit_" + setName_) +\
                            estimateH2Simple.only_nodes_with_namespace("herit_" + setName_)
    setName_ = 'TwentyTwoSets'
    # twentyTwoSetsChromosomes = setAnalysisEnvironmentPipeline + \
    #                         calculateZZt.only_nodes_with_namespace("herit_" + setName_) +\
    #                         estimateH2Simple.only_nodes_with_namespace("herit_" + setName_)


    return {
        "pipelineFullChrSer": pipelineFullChrSer,
        "twoSetsChromosomes": twoSetsChromosomes,
        # "twentyTwoSetsChromosomes": twentyTwoSetsChromosomes,
        "__default__": setAnalysisEnvironmentPipeline + calculateZZt +  estimateH2Simple 
    }
