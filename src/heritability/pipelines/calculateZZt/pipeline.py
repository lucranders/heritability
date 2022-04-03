"""
This is a boilerplate pipeline 'calculateZZt'
generated using Kedro 0.17.6
"""

from kedro.pipeline import Pipeline, node
from .nodes import calculateGCTA, correctGRM

def create_pipeline(**kwargs):
    return Pipeline([
        node(
            calculateGCTA,
            ["params:nameFile",'params:listChrs','params:nameMatrix',"params:tempFolderPath","params:gcta","params:numThreads","bedFiles"],
            outputs='calculatedZZt',
            name="calculateZZt",
        ),
        node(
            correctGRM,
            ["params:tempFolderPath",'params:nameMatrix',"calculatedZZt"],
            outputs='correctedZZt',
            name="ZZtCorrection",
        ),
    ])
