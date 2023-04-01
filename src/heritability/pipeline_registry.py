"""Project pipelines."""
from typing import Dict

from kedro.pipeline import Pipeline

from heritability.pipelines import _01_set_analysis_environment as sae
from heritability.pipelines import _02_calculateZZt as cZZt
from heritability.pipelines import _03_estimate_h2 as eH2S


def register_pipelines() -> Dict[str, Pipeline]:
    """Register the project's pipelines.

    Returns:
        A mapping from a pipeline name to a ``Pipeline`` object.

    """
    setAnalysisEnvironmentPipeline = sae.create_pipeline()
    calculateZZt = cZZt.create_pipeline()
    estimate_heritability = eH2S.create_pipeline()

    full_pipeline = setAnalysisEnvironmentPipeline +\
                            calculateZZt +\
                            estimate_heritability


    return {
        "__default__": full_pipeline
    }
