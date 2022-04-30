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
            ["params:nameFile",'params:listChrs','params:nameMatrix',"params:pathTempFiles","params:numThreads","bed_Files"],
            outputs='calculated_ZZt',
            name="calculate_ZZt",
        ),
        node(
            correctGRM,
            ["params:pathTempFiles",'params:nameMatrix',"calculated_ZZt"],
            outputs='corrected_ZZt',
            name="ZZt_Correction",
        ),
    ])
