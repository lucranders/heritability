"""Project pipelines."""
from typing import Dict

from kedro.pipeline import Pipeline

from heritability.pipelines import data_engineering as de
from heritability.pipelines import data_science as ds
from heritability.pipelines import setAnalysisEnvironment as sae


def register_pipelines() -> Dict[str, Pipeline]:
    """Register the project's pipelines.

    Returns:
        A mapping from a pipeline name to a ``Pipeline`` object.

    """
    data_engineering_pipeline = de.create_pipeline()
    data_science_pipeline = ds.create_pipeline()
    setAnalysisEnvironmentPipeline = sae.create_pipeline()

    return {
        # "de": data_engineering_pipeline,
        # "ds": data_science_pipeline,
        "sae": setAnalysisEnvironmentPipeline,
        # "__default__": setAnalysisEnvironmentPipeline + data_engineering_pipeline + data_science_pipeline,
        "__default__": setAnalysisEnvironmentPipeline
    }
