"""Project pipelines."""
from typing import Dict

from kedro.pipeline import Pipeline

from heritability.pipelines import setAnalysisEnvironment as sae
from heritability.pipelines import calculateZZt as cZZt
# from heritability.pipelines import calculateZZtV2 as cZZtV2
from heritability.pipelines import estimateH2Simple as eH2S


def register_pipelines() -> Dict[str, Pipeline]:
    """Register the project's pipelines.

    Returns:
        A mapping from a pipeline name to a ``Pipeline`` object.

    """
    setAnalysisEnvironmentPipeline = sae.create_pipeline()
    calculateZZt = cZZt.create_pipeline()
    # calculateZZtV2 = cZZtV2.create_pipeline()
    estimateH2Simple = eH2S.create_pipeline()

    return {
        "sae": setAnalysisEnvironmentPipeline,
        "cZZt":calculateZZt,
        # "cZZtV2":calculateZZtV2,
        "eH2S":estimateH2Simple,
        "__default__": setAnalysisEnvironmentPipeline + calculateZZt + estimateH2Simple
    }
